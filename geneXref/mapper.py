"""GeneMapper — map gene and transcript identifiers via a SQLite backend."""

from __future__ import annotations

import sqlite3
import warnings
from pathlib import Path

import pandas as pd


# ID types that live in the gene-level views/tables
GENE_ID_TYPES = {
    "ensembl_gene_id",
    "gene_name",
    "ncbi_gene_id",
    "hgnc_id",
}

# ID types that live in the transcript-level views/tables
TRANSCRIPT_ID_TYPES = {
    "ensembl_transcript_id",
    "refseq_mrna_id",
    "ucsc_transcript_id",
}

ALL_ID_TYPES = sorted(GENE_ID_TYPES | TRANSCRIPT_ID_TYPES)

# Which view to query for a given (input, output) pair.
# When both IDs are gene-level, use ensg2all.
# When both are transcript-level, use enst2all.
# Cross-level queries join the two views.
_GENE_VIEW = "ensg2all"
_TRANSCRIPT_VIEW = "enst2all"


def _default_db_path() -> Path:
    """Locate the database: check ~/.geneXref/ first, then the package data dir."""
    user_db = Path.home() / ".geneXref" / "geneXref.db"
    if user_db.exists():
        return user_db
    bundled = Path(__file__).parent / "data" / "geneXref.db"
    if bundled.exists():
        return bundled
    raise FileNotFoundError(
        "No geneXref database found. Either:\n"
        "  - Run: python -m geneXref.download  (to fetch a pre-built database)\n"
        "  - Run: python build_db.py           (to build from Ensembl)\n"
        "  - Pass db_path= to GeneMapper()"
    )


class GeneMapper:
    """Map gene and transcript identifiers across Ensembl, NCBI, HGNC, RefSeq, and UCSC.

    Parameters
    ----------
    db_path : str or Path, optional
        Path to a SQLite database built by ``build_db.py``.  If *None*,
        looks for ``~/.geneXref/geneXref.db`` (downloaded), then falls back
        to any database bundled with the package.
    """

    def __init__(self, db_path: str | Path | None = None):
        if db_path is None:
            db_path = _default_db_path()
        self._db_path = Path(db_path)
        if not self._db_path.exists():
            raise FileNotFoundError(
                f"Database not found: {self._db_path}. "
                "Run 'python -m geneXref.download' or build_db.py first."
            )
        self._conn = sqlite3.connect(str(self._db_path))
        self._conn.execute("PRAGMA foreign_keys = ON")

    def close(self):
        """Close the underlying database connection."""
        self._conn.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    # ------------------------------------------------------------------ #
    # Public API
    # ------------------------------------------------------------------ #

    @staticmethod
    def list_id_types() -> list[str]:
        """Return the list of supported identifier types."""
        return list(ALL_ID_TYPES)

    def map(
        self,
        ids: list[str],
        input_id: str,
        output_id: str,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Map a list of identifiers from one type to another.

        Parameters
        ----------
        ids : list[str]
            Identifiers to map.
        input_id : str
            Source identifier type (see :meth:`list_id_types`).
        output_id : str
            Target identifier type (see :meth:`list_id_types`).

        Returns
        -------
        mapped : DataFrame
            Successfully mapped identifiers with columns ``[input_id, output_id]``.
        unmapped : DataFrame
            Identifiers that could not be mapped, with columns
            ``["identifier", "reason"]``.  Reasons: ``"not_found"``,
            ``"duplicate_input"``, ``"duplicate_output"``.
        """
        self._validate_id_type(input_id)
        self._validate_id_type(output_id)
        if input_id == output_id:
            raise ValueError("input_id and output_id must be different")

        # Strip Ensembl version suffixes for lookup, but keep originals
        original_ids = list(ids)
        lookup_ids = [self._strip_version(i) for i in original_ids]
        lookup_to_original = dict(zip(lookup_ids, original_ids))

        # Query the database
        view, columns = self._resolve_query(input_id, output_id)
        input_col, output_col = columns
        raw = self._query_mapping(view, input_col, output_col, lookup_ids)

        # Build result DataFrames
        mapped_rows: list[tuple[str, str]] = []
        unmapped_rows: list[tuple[str, str]] = []

        # Detect duplicates
        input_counts: dict[str, list[str]] = {}
        output_counts: dict[str, list[str]] = {}
        for inp, out in raw:
            input_counts.setdefault(inp, []).append(out)
            output_counts.setdefault(out, []).append(inp)

        dup_inputs = {k for k, v in input_counts.items() if len(v) > 1}
        dup_outputs = {k for k, v in output_counts.items() if len(v) > 1}

        seen_inputs: set[str] = set()
        for inp, out in raw:
            if inp in dup_inputs:
                if inp not in seen_inputs:
                    unmapped_rows.append((lookup_to_original.get(inp, inp), "duplicate_input"))
                    seen_inputs.add(inp)
            elif out in dup_outputs:
                pass  # handled below
            else:
                mapped_rows.append((lookup_to_original.get(inp, inp), out))

        # Report duplicate outputs (multiple inputs -> same output)
        for out, inputs in output_counts.items():
            if len(inputs) > 1:
                for inp in inputs:
                    if inp not in dup_inputs:
                        unmapped_rows.append((lookup_to_original.get(inp, inp), "duplicate_output"))

        # Report not-found
        found_inputs = {row[0] for row in raw}
        for lid, oid in zip(lookup_ids, original_ids):
            if lid not in found_inputs:
                unmapped_rows.append((oid, "not_found"))

        mapped = pd.DataFrame(mapped_rows, columns=[input_id, output_id])
        unmapped = pd.DataFrame(
            unmapped_rows, columns=["identifier", "reason"]
        ) if unmapped_rows else pd.DataFrame(columns=["identifier", "reason"])

        if len(unmapped) > 0:
            n = len(unmapped)
            warnings.warn(f"{n} identifier(s) could not be mapped", UserWarning, stacklevel=2)

        return mapped, unmapped

    # ------------------------------------------------------------------ #
    # Internals
    # ------------------------------------------------------------------ #

    @staticmethod
    def _validate_id_type(id_type: str):
        if id_type not in ALL_ID_TYPES:
            raise ValueError(
                f"Unknown id_type {id_type!r}. "
                f"Valid types: {ALL_ID_TYPES}"
            )

    @staticmethod
    def _strip_version(identifier: str) -> str:
        """Strip Ensembl version suffixes like ENSG00000141510.14."""
        if identifier.startswith(("ENSG", "ENST")) and "." in identifier:
            return identifier.split(".")[0]
        return identifier

    def _resolve_query(
        self, input_id: str, output_id: str
    ) -> tuple[str, tuple[str, str]]:
        """Determine which view/join to query and the column names to select.

        Returns (view_or_subquery, (input_column, output_column)).
        """
        in_gene = input_id in GENE_ID_TYPES
        out_gene = output_id in GENE_ID_TYPES

        if in_gene and out_gene:
            return _GENE_VIEW, (input_id, output_id)

        if not in_gene and not out_gene:
            return _TRANSCRIPT_VIEW, (input_id, output_id)

        # Cross-level: join gene and transcript views on ensembl_gene_id
        subquery = (
            f"(SELECT g.*, t.ensembl_transcript_id, t.transcript_biotype, "
            f"t.refseq_mrna_id, t.ucsc_transcript_id "
            f"FROM {_GENE_VIEW} g "
            f"JOIN {_TRANSCRIPT_VIEW} t "
            f"ON g.ensembl_gene_id = t.ensembl_gene_id)"
        )
        return subquery, (input_id, output_id)

    def _query_mapping(
        self,
        source: str,
        input_col: str,
        output_col: str,
        ids: list[str],
    ) -> list[tuple[str, str]]:
        """Query the database for (input, output) pairs matching the given IDs."""
        if not ids:
            return []

        # SQLite has a variable limit; batch in chunks of 500
        results: list[tuple[str, str]] = []
        for i in range(0, len(ids), 500):
            chunk = ids[i : i + 500]
            placeholders = ",".join("?" * len(chunk))
            sql = (
                f"SELECT {input_col}, {output_col} FROM {source} "
                f"WHERE {input_col} IN ({placeholders}) "
                f"AND {output_col} IS NOT NULL"
            )
            results.extend(self._conn.execute(sql, chunk).fetchall())

        return results
