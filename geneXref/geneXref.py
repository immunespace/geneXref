import warnings
import pandas as pd


class geneXref:
    """Maps gene identifiers across databases using a pre-built TSV database."""

    def __init__(self, db_path: str):
        """Load a geneXref database from a TSV file.

        Parameters
        ----------
        db_path : str
            Path to the geneXref database TSV file.
        """
        self._db = pd.read_csv(db_path, sep="\t", dtype=str)

    def list_id_types(self) -> list[str]:
        """Return the list of identifier types available in the database."""
        return list(self._db.columns)

    def map(
        self,
        ids: list[str],
        input_id: str,
        output_ids: list[str],
        remove_unmapped: bool = False,
    ) -> pd.DataFrame:
        """Map a list of identifiers to one or more other identifier types.

        A mapping is set to NaN and a warning is issued when:
        - No row matches the input identifier.
        - More than one row matches the input identifier.
        - The output value found is shared by a different row (many-to-many).

        Parameters
        ----------
        ids : list[str]
            Input identifiers to map.
        input_id : str
            Column name of the input identifier type.
        output_ids : list[str]
            Column names of the desired output identifier types.
        remove_unmapped : bool, optional
            If True, drop rows where any output identifier could not be
            resolved (unmapped or ambiguous).  Defaults to False.

        Returns
        -------
        pd.DataFrame
            DataFrame with columns [input_id] + output_ids.
        """
        if input_id not in self._db.columns:
            raise ValueError(
                f"'{input_id}' is not a recognised identifier type. "
                f"Call list_id_types() to see available types."
            )
        missing = [c for c in output_ids if c not in self._db.columns]
        if missing:
            raise ValueError(
                f"The following output identifier types are not recognised: "
                f"{missing}. Call list_id_types() to see available types."
            )

        # Pre-compute how many times each value appears in each output column
        # so we can flag output values that are shared across multiple rows.
        output_value_counts = {
            col: self._db[col].value_counts() for col in output_ids
        }

        result_rows = []
        for query_id in ids:
            hits = self._db[self._db[input_id] == query_id]

            if len(hits) == 0:
                warnings.warn(
                    f"No mapping found for {input_id}='{query_id}'.",
                    stacklevel=2,
                )
                row = {input_id: query_id, **{col: pd.NA for col in output_ids}}

            elif len(hits) > 1:
                warnings.warn(
                    f"Multiple rows found for {input_id}='{query_id}'; "
                    f"all output values set to NaN.",
                    stacklevel=2,
                )
                row = {input_id: query_id, **{col: pd.NA for col in output_ids}}

            else:
                row = {input_id: query_id}
                for col in output_ids:
                    val = hits.iloc[0][col]
                    if pd.isna(val):
                        row[col] = pd.NA
                    elif output_value_counts[col].get(val, 0) > 1:
                        warnings.warn(
                            f"Output value '{val}' in '{col}' is linked to multiple "
                            f"entries; setting to NaN for {input_id}='{query_id}'.",
                            stacklevel=2,
                        )
                        row[col] = pd.NA
                    else:
                        row[col] = val

            result_rows.append(row)

        result = pd.DataFrame(result_rows, columns=[input_id] + output_ids)

        if remove_unmapped:
            result = result.dropna(subset=output_ids).reset_index(drop=True)

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
        df["ensembl_transcript_id"] = df["mane_select"].str.split("|").str[0]
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
