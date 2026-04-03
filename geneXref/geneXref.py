import re
import warnings
from importlib import resources

import pandas as pd

# Matches versioned Ensembl identifiers, e.g. ENSG00000139618.14
_ENSEMBL_VERSION_RE = re.compile(r"^(ENS[A-Z]*\d+)\.\d+$")


class geneXref:
    """Maps gene identifiers across databases using a pre-built TSV database."""

    def __init__(self, db_path: str | None = None):
        """Load a geneXref database from a TSV file.

        Parameters
        ----------
        db_path : str, optional
            Path to a geneXref database TSV file.  When *None* (the
            default), the pre-built database bundled with the package is
            used.
        """
        if db_path is None:
            db_path = str(
                resources.files("geneXref").joinpath("data", "db-20260403.tsv")
            )
        self._db = pd.read_csv(db_path, sep="\t", dtype=str)

    def list_id_types(self) -> list[str]:
        """Return the list of identifier types available in the database."""
        return list(self._db.columns)

    def map(
        self,
        ids: list[str],
        input_id: str,
        output_id: str,
        remove_unmapped: bool = False,
    ) -> pd.DataFrame:
        """Map a list of identifiers to another identifier type.

        Parameters
        ----------
        ids : list[str]
            Input identifiers to map.
        input_id : str
            Column name of the input identifier type.
        output_id : str
            Column name of the desired output identifier type.
        remove_unmapped : bool, optional
            If True, drop rows where the output identifier could not be
            resolved.  Defaults to False.

        Returns
        -------
        pd.DataFrame
            DataFrame with columns [input_id, output_id].
        """
        if input_id not in self._db.columns:
            raise ValueError(
                f"'{input_id}' is not a recognised identifier type. "
                f"Call list_id_types() to see available types."
            )
        if output_id not in self._db.columns:
            raise ValueError(
                f"'{output_id}' is not a recognised identifier type. "
                f"Call list_id_types() to see available types."
            )

        # Strip Ensembl version suffixes from input IDs
        lookup_ids = pd.Series(ids).str.replace(
            _ENSEMBL_VERSION_RE, r"\1", regex=True
        )

        # Extract only the two relevant columns; drop missing output values
        sub = self._db[[input_id, output_id]].dropna(subset=[output_id])

        # Identify ambiguous output values across the full database
        dup_outputs = sub[output_id].duplicated(keep=False)
        sub = sub[~dup_outputs]

        # Filter to rows matching the input IDs
        sub = sub[sub[input_id].isin(lookup_ids)]

        # Remove duplicate input rows (ambiguous input mappings)
        dup_inputs = sub[input_id].duplicated(keep=False)
        sub = sub[~dup_inputs]

        # Build the result: start from the original IDs, merge with matches
        query = pd.DataFrame({input_id: ids, "_lookup": lookup_ids})
        result = query.merge(sub, left_on="_lookup", right_on=input_id,
                             how="left", suffixes=("", "_db"))
        result = result[[input_id, output_id]]

        n_unmapped = result[output_id].isna().sum()
        if n_unmapped > 0:
            warnings.warn(
                f"{n_unmapped} of {len(ids)} identifiers could not be mapped.",
                stacklevel=2,
            )

        if remove_unmapped:
            result = result.dropna(subset=[output_id]).reset_index(drop=True)

        return result

    @staticmethod
    def rebuild_database(hgnc_path: str, output_path: str) -> None:
        """Build a geneXref database from an HGNC complete-set export.

        Columns retained from the HGNC file and how they are processed:

        - hgnc_id              : numeric portion stripped of the 'HGNC:' prefix
        - symbol               : renamed to gene_symbol
        - entrez_id            : renamed to ncbi_gene_id
        - ensembl_gene_id      : kept as-is
        - ucsc_id              : kept as-is
        - refseq_accession     : kept as-is
        - mane_select          : pipe-delimited; Ensembl transcript ID (first
                                 field) extracted and saved as ensembl_transcript_id
        - uniprot_ids          : pipe-delimited; first entry kept and renamed
                                 to uniprot_id

        Parameters
        ----------
        hgnc_path : str
            Path to the downloaded HGNC complete set TSV file.
        output_path : str
            Destination path for the geneXref database TSV.
        """
        KEEP_COLS = [
            "hgnc_id",
            "symbol",
            "entrez_id",
            "ensembl_gene_id",
            "ucsc_id",
            "refseq_accession",
            "mane_select",
            "uniprot_ids",
        ]

        df = pd.read_csv(hgnc_path, sep="\t", dtype=str, usecols=KEEP_COLS)

        # hgnc_id: strip 'HGNC:' prefix, keep numeric portion
        df["hgnc_id"] = df["hgnc_id"].str.replace(r"^HGNC:", "", regex=True)

        # mane_select: pipe-delimited "ENST…|NM_…"; keep Ensembl transcript ID
        # and strip the version suffix so the DB stores bare IDs consistently.
        df["ensembl_transcript_id"] = (
            df["mane_select"].str.split("|").str[0]
            .str.replace(r"\.\d+$", "", regex=True)
        )
        df = df.drop(columns=["mane_select"])

        # uniprot_ids: keep only the first pipe-delimited entry
        df["uniprot_ids"] = df["uniprot_ids"].str.split("|").str[0]

        df = df.rename(columns={
            "symbol": "gene_symbol",
            "entrez_id": "ncbi_gene_id",
            "uniprot_ids": "uniprot_id",
        })

        df = df[[
            "hgnc_id",
            "gene_symbol",
            "ncbi_gene_id",
            "ensembl_gene_id",
            "ucsc_id",
            "refseq_accession",
            "ensembl_transcript_id",
            "uniprot_id",
        ]]

        df.to_csv(output_path, sep="\t", index=False)
