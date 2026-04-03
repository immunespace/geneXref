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
            DataFrame with columns [input_id] + output_ids.  Rows that could
            not be mapped unambiguously contain NaN in the output columns.
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
                    f"Multiple mappings found for {input_id}='{query_id}'; "
                    f"output set to NaN.",
                    stacklevel=2,
                )
                row = {input_id: query_id, **{col: pd.NA for col in output_ids}}

            else:
                row = {input_id: query_id}
                for col in output_ids:
                    row[col] = hits.iloc[0][col]

            result_rows.append(row)

        result = pd.DataFrame(result_rows, columns=[input_id] + output_ids)

        if remove_unmapped:
            result = result.dropna(subset=output_ids).reset_index(drop=True)

        return result

    @staticmethod
    def rebuild_database(hgnc_path: str, output_path: str) -> None:
        """Build a geneXref database from an HGNC complete-set export.

        Parameters
        ----------
        hgnc_path : str
            Path to the downloaded HGNC complete set TSV file.
        output_path : str
            Destination path for the geneXref database TSV.
        """
        raise NotImplementedError("rebuild_database is not yet implemented.")
