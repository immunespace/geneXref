import warnings

import pytest
import pandas as pd

from geneXref import geneXref


# ===========================================================================
# rebuild_database
# ===========================================================================

class TestRebuildDatabase:

    def test_output_columns(self, hgnc_path, tmp_path):
        out = str(tmp_path / "db.tsv")
        geneXref.rebuild_database(hgnc_path, out)
        df = pd.read_csv(out, sep="\t", dtype=str)
        assert list(df.columns) == [
            "hgnc_id", "gene_symbol", "ncbi_gene_id", "ensembl_gene_id",
            "ucsc_id", "refseq_accession", "ensembl_transcript_id", "uniprot_id",
        ]

    def test_hgnc_id_prefix_stripped(self, hgnc_path, tmp_path):
        out = str(tmp_path / "db.tsv")
        geneXref.rebuild_database(hgnc_path, out)
        df = pd.read_csv(out, sep="\t", dtype=str)
        assert not df["hgnc_id"].str.startswith("HGNC:").any()
        assert "1944" in df["hgnc_id"].values

    def test_column_renames(self, hgnc_path, tmp_path):
        out = str(tmp_path / "db.tsv")
        geneXref.rebuild_database(hgnc_path, out)
        df = pd.read_csv(out, sep="\t", dtype=str)
        assert "gene_symbol" in df.columns
        assert "ncbi_gene_id" in df.columns
        assert "uniprot_id" in df.columns
        # original names must be absent
        assert "symbol" not in df.columns
        assert "entrez_id" not in df.columns
        assert "uniprot_ids" not in df.columns
        assert "mane_select" not in df.columns

    def test_mane_select_extracts_ensembl_transcript(self, hgnc_path, tmp_path):
        out = str(tmp_path / "db.tsv")
        geneXref.rebuild_database(hgnc_path, out)
        df = pd.read_csv(out, sep="\t", dtype=str)
        chn2 = df[df["gene_symbol"] == "CHN2"].iloc[0]
        assert chn2["ensembl_transcript_id"] == "ENST00000222792"

    def test_mane_select_refseq_portion_discarded(self, hgnc_path, tmp_path):
        out = str(tmp_path / "db.tsv")
        geneXref.rebuild_database(hgnc_path, out)
        df = pd.read_csv(out, sep="\t", dtype=str)
        # No NM_ accession should appear in the ensembl_transcript_id column
        has_nm = df["ensembl_transcript_id"].dropna().str.startswith("NM_")
        assert not has_nm.any()

    def test_uniprot_first_entry_only(self, hgnc_path, tmp_path):
        out = str(tmp_path / "db.tsv")
        geneXref.rebuild_database(hgnc_path, out)
        df = pd.read_csv(out, sep="\t", dtype=str)
        map2k1 = df[df["gene_symbol"] == "MAP2K1"].iloc[0]
        # Original had "Q02750|A0A0A0MRZ7"; only first should be kept
        assert map2k1["uniprot_id"] == "Q02750"

    def test_missing_optional_fields_are_nan(self, hgnc_path, tmp_path):
        out = str(tmp_path / "db.tsv")
        geneXref.rebuild_database(hgnc_path, out)
        df = pd.read_csv(out, sep="\t", dtype=str)
        ska2p1 = df[df["gene_symbol"] == "SKA2P1"].iloc[0]
        assert pd.isna(ska2p1["ensembl_transcript_id"])
        assert pd.isna(ska2p1["uniprot_id"])
        assert pd.isna(ska2p1["ucsc_id"])

    def test_row_count_preserved(self, hgnc_path, tmp_path):
        out = str(tmp_path / "db.tsv")
        geneXref.rebuild_database(hgnc_path, out)
        df = pd.read_csv(out, sep="\t", dtype=str)
        assert len(df) == 5


# ===========================================================================
# list_id_types
# ===========================================================================

class TestListIdTypes:

    def test_returns_all_columns(self, db_path):
        gx = geneXref(db_path)
        assert gx.list_id_types() == [
            "hgnc_id", "gene_symbol", "ncbi_gene_id", "ensembl_gene_id",
            "ucsc_id", "refseq_accession", "ensembl_transcript_id", "uniprot_id",
        ]


# ===========================================================================
# map
# ===========================================================================

class TestMap:

    # --- successful mapping -------------------------------------------------

    def test_single_id_maps_correctly(self, db_path):
        gx = geneXref(db_path)
        result = gx.map(["ENSG00000106069"],
                        input_id="ensembl_gene_id",
                        output_id="gene_symbol")
        assert result.loc[0, "gene_symbol"] == "CHN2"

    def test_multiple_ids_map_correctly(self, db_path):
        gx = geneXref(db_path)
        result = gx.map(["ENSG00000106069", "ENSG00000187753"],
                        input_id="ensembl_gene_id",
                        output_id="gene_symbol")
        assert list(result["gene_symbol"]) == ["CHN2", "C9orf153"]

    def test_output_dataframe_has_correct_columns(self, db_path):
        gx = geneXref(db_path)
        result = gx.map(["ENSG00000106069"],
                        input_id="ensembl_gene_id",
                        output_id="gene_symbol")
        assert list(result.columns) == ["ensembl_gene_id", "gene_symbol"]

    # --- unmapped IDs -------------------------------------------------------

    def test_unmapped_id_produces_nan(self, db_path):
        gx = geneXref(db_path)
        with pytest.warns(UserWarning, match="could not be mapped"):
            result = gx.map(["ENSG_DOES_NOT_EXIST"],
                            input_id="ensembl_gene_id",
                            output_id="gene_symbol")
        assert pd.isna(result.loc[0, "gene_symbol"])

    def test_unmapped_id_warns_with_count(self, db_path):
        gx = geneXref(db_path)
        with pytest.warns(UserWarning, match="1 of 1 identifiers could not be mapped"):
            gx.map(["ENSG_DOES_NOT_EXIST"],
                   input_id="ensembl_gene_id",
                   output_id="gene_symbol")

    # --- duplicate input IDs ------------------------------------------------

    def test_duplicate_input_id_produces_nan(self, db_path, tmp_path):
        # Write a db where the same ensembl_gene_id appears twice
        import pandas as pd
        df = pd.read_csv(db_path, sep="\t", dtype=str)
        dup = df[df["ensembl_gene_id"] == "ENSG00000106069"].copy()
        dup["gene_symbol"] = "CHN2_DUP"
        df = pd.concat([df, dup], ignore_index=True)
        dup_path = str(tmp_path / "dup_db.tsv")
        df.to_csv(dup_path, sep="\t", index=False)

        gx = geneXref(dup_path)
        with pytest.warns(UserWarning, match="could not be mapped"):
            result = gx.map(["ENSG00000106069"],
                            input_id="ensembl_gene_id",
                            output_id="gene_symbol")
        assert pd.isna(result.loc[0, "gene_symbol"])

    # --- many-to-many (shared output value) ---------------------------------

    def test_shared_output_value_produces_nan(self, db_path):
        # uniprot_id "P52757" is shared by CHN2 and DUPGENE in the test db
        gx = geneXref(db_path)
        with pytest.warns(UserWarning, match="could not be mapped"):
            result = gx.map(["ENSG00000106069"],
                            input_id="ensembl_gene_id",
                            output_id="uniprot_id")
        assert pd.isna(result.loc[0, "uniprot_id"])

    def test_unambiguous_output_is_not_warned(self, db_path):
        gx = geneXref(db_path)
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            # gene_symbol is unique — should raise no warning
            gx.map(["ENSG00000106069"],
                   input_id="ensembl_gene_id",
                   output_id="gene_symbol")

    # --- remove_unmapped ----------------------------------------------------

    def test_remove_unmapped_drops_nan_rows(self, db_path):
        gx = geneXref(db_path)
        with pytest.warns(UserWarning):
            result = gx.map(["ENSG00000106069", "ENSG_FAKE"],
                            input_id="ensembl_gene_id",
                            output_id="gene_symbol",
                            remove_unmapped=True)
        assert len(result) == 1
        assert result.loc[0, "gene_symbol"] == "CHN2"

    def test_remove_unmapped_false_keeps_nan_rows(self, db_path):
        gx = geneXref(db_path)
        with pytest.warns(UserWarning):
            result = gx.map(["ENSG00000106069", "ENSG_FAKE"],
                            input_id="ensembl_gene_id",
                            output_id="gene_symbol",
                            remove_unmapped=False)
        assert len(result) == 2

    # --- validation ---------------------------------------------------------

    def test_invalid_input_id_raises(self, db_path):
        gx = geneXref(db_path)
        with pytest.raises(ValueError, match="not a recognised identifier type"):
            gx.map(["X"], input_id="nonexistent_col", output_id="gene_symbol")

    def test_invalid_output_id_raises(self, db_path):
        gx = geneXref(db_path)
        with pytest.raises(ValueError, match="not a recognised identifier type"):
            gx.map(["ENSG00000106069"],
                   input_id="ensembl_gene_id",
                   output_id="nonexistent_col")

    # --- Ensembl version stripping ------------------------------------------

    def test_versioned_ensembl_gene_id_maps(self, db_path):
        gx = geneXref(db_path)
        result = gx.map(["ENSG00000106069.7"],
                        input_id="ensembl_gene_id",
                        output_id="gene_symbol")
        assert result.loc[0, "gene_symbol"] == "CHN2"

    def test_versioned_id_preserved_in_output(self, db_path):
        # The original versioned value should appear in the input_id column
        gx = geneXref(db_path)
        result = gx.map(["ENSG00000106069.7"],
                        input_id="ensembl_gene_id",
                        output_id="gene_symbol")
        assert result.loc[0, "ensembl_gene_id"] == "ENSG00000106069.7"

    def test_versioned_ensembl_transcript_id_maps(self, db_path):
        # DB stores unversioned IDs; versioned query should still resolve.
        gx = geneXref(db_path)
        result = gx.map(["ENST00000222792.11"],
                        input_id="ensembl_transcript_id",
                        output_id="gene_symbol")
        assert result.loc[0, "gene_symbol"] == "CHN2"

    def test_non_ensembl_id_unaffected_by_version_stripping(self, db_path):
        # Gene symbols containing dots should not be altered
        gx = geneXref(db_path)
        with pytest.warns(UserWarning, match="could not be mapped"):
            result = gx.map(["FAKEGENE.1"],
                            input_id="gene_symbol",
                            output_id="ensembl_gene_id")
        assert pd.isna(result.loc[0, "ensembl_gene_id"])
