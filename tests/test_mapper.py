"""Tests for GeneMapper."""

import pytest
from geneXref import GeneMapper


class TestListIdTypes:
    def test_returns_all_types(self):
        types = GeneMapper.list_id_types()
        assert "ensembl_gene_id" in types
        assert "gene_name" in types
        assert "ncbi_gene_id" in types
        assert "hgnc_id" in types
        assert "ensembl_transcript_id" in types
        assert "refseq_mrna_id" in types
        assert "ucsc_transcript_id" in types

    def test_returns_sorted(self):
        types = GeneMapper.list_id_types()
        assert types == sorted(types)


class TestMapGeneLevel:
    def test_ensg_to_ncbi(self, test_db):
        with GeneMapper(test_db) as gm:
            mapped, unmapped = gm.map(
                ["ENSG00000141510", "ENSG00000012048"],
                input_id="ensembl_gene_id",
                output_id="ncbi_gene_id",
            )
        assert len(mapped) == 2
        assert len(unmapped) == 0
        row = mapped[mapped["ensembl_gene_id"] == "ENSG00000141510"]
        assert row["ncbi_gene_id"].iloc[0] == "7157"

    def test_ncbi_to_ensg(self, test_db):
        with GeneMapper(test_db) as gm:
            mapped, unmapped = gm.map(
                ["7157", "672"],
                input_id="ncbi_gene_id",
                output_id="ensembl_gene_id",
            )
        assert len(mapped) == 2
        row = mapped[mapped["ncbi_gene_id"] == "672"]
        assert row["ensembl_gene_id"].iloc[0] == "ENSG00000012048"

    def test_gene_name_to_ensg(self, test_db):
        with GeneMapper(test_db) as gm:
            mapped, unmapped = gm.map(
                ["TP53", "BRCA1"],
                input_id="gene_name",
                output_id="ensembl_gene_id",
            )
        assert len(mapped) == 2
        assert len(unmapped) == 0

    def test_ensg_to_hgnc(self, test_db):
        with GeneMapper(test_db) as gm:
            mapped, unmapped = gm.map(
                ["ENSG00000141510"],
                input_id="ensembl_gene_id",
                output_id="hgnc_id",
            )
        assert mapped["hgnc_id"].iloc[0] == "HGNC:11998"

    def test_ncbi_to_hgnc_via_ensg(self, test_db):
        """Indirect mapping: ncbi -> hgnc goes through ensg2all view."""
        with GeneMapper(test_db) as gm:
            mapped, unmapped = gm.map(
                ["7157"],
                input_id="ncbi_gene_id",
                output_id="hgnc_id",
            )
        assert len(mapped) == 1
        assert mapped["hgnc_id"].iloc[0] == "HGNC:11998"


class TestMapTranscriptLevel:
    def test_enst_to_refseq(self, test_db):
        with GeneMapper(test_db) as gm:
            mapped, unmapped = gm.map(
                ["ENST00000269305"],
                input_id="ensembl_transcript_id",
                output_id="refseq_mrna_id",
            )
        assert mapped["refseq_mrna_id"].iloc[0] == "NM_000546"

    def test_refseq_to_enst(self, test_db):
        with GeneMapper(test_db) as gm:
            mapped, unmapped = gm.map(
                ["NM_007294"],
                input_id="refseq_mrna_id",
                output_id="ensembl_transcript_id",
            )
        assert mapped["ensembl_transcript_id"].iloc[0] == "ENST00000357654"

    def test_enst_to_ucsc(self, test_db):
        with GeneMapper(test_db) as gm:
            mapped, unmapped = gm.map(
                ["ENST00000275493"],
                input_id="ensembl_transcript_id",
                output_id="ucsc_transcript_id",
            )
        assert mapped["ucsc_transcript_id"].iloc[0] == "uc003tpa.3"


class TestMapCrossLevel:
    def test_gene_name_to_refseq(self, test_db):
        """Cross-level: gene name -> refseq transcript."""
        with GeneMapper(test_db) as gm:
            mapped, unmapped = gm.map(
                ["TP53"],
                input_id="gene_name",
                output_id="refseq_mrna_id",
            )
        # TP53 has ENST00000269305 -> NM_000546, and ENST00000413465 with no refseq
        assert "NM_000546" in mapped["refseq_mrna_id"].values

    def test_refseq_to_gene_name(self, test_db):
        """Cross-level: refseq -> gene name."""
        with GeneMapper(test_db) as gm:
            mapped, unmapped = gm.map(
                ["NM_000546"],
                input_id="refseq_mrna_id",
                output_id="gene_name",
            )
        assert mapped["gene_name"].iloc[0] == "TP53"


class TestUnmapped:
    def test_not_found(self, test_db):
        with GeneMapper(test_db) as gm:
            mapped, unmapped = gm.map(
                ["ENSG_FAKE_ID"],
                input_id="ensembl_gene_id",
                output_id="ncbi_gene_id",
            )
        assert len(mapped) == 0
        assert len(unmapped) == 1
        assert unmapped["reason"].iloc[0] == "not_found"

    def test_mixed_found_and_not_found(self, test_db):
        with GeneMapper(test_db) as gm:
            mapped, unmapped = gm.map(
                ["ENSG00000141510", "ENSG_DOES_NOT_EXIST"],
                input_id="ensembl_gene_id",
                output_id="ncbi_gene_id",
            )
        assert len(mapped) == 1
        assert len(unmapped) == 1
        assert unmapped["identifier"].iloc[0] == "ENSG_DOES_NOT_EXIST"

    def test_no_mapping_for_existing_gene(self, test_db):
        """LNCRNA1 exists in ensg but has no NCBI mapping."""
        with GeneMapper(test_db) as gm:
            mapped, unmapped = gm.map(
                ["ENSG00000099999"],
                input_id="ensembl_gene_id",
                output_id="ncbi_gene_id",
            )
        assert len(mapped) == 0
        assert len(unmapped) == 1
        assert unmapped["reason"].iloc[0] == "not_found"

    def test_warning_issued(self, test_db):
        with GeneMapper(test_db) as gm:
            with pytest.warns(UserWarning, match="could not be mapped"):
                gm.map(
                    ["FAKE"],
                    input_id="ensembl_gene_id",
                    output_id="ncbi_gene_id",
                )


class TestVersionStripping:
    def test_ensg_version_stripped(self, test_db):
        with GeneMapper(test_db) as gm:
            mapped, unmapped = gm.map(
                ["ENSG00000141510.14"],
                input_id="ensembl_gene_id",
                output_id="ncbi_gene_id",
            )
        assert len(mapped) == 1
        # Original versioned ID preserved in output
        assert mapped["ensembl_gene_id"].iloc[0] == "ENSG00000141510.14"

    def test_enst_version_stripped(self, test_db):
        with GeneMapper(test_db) as gm:
            mapped, unmapped = gm.map(
                ["ENST00000269305.8"],
                input_id="ensembl_transcript_id",
                output_id="refseq_mrna_id",
            )
        assert len(mapped) == 1
        assert mapped["ensembl_transcript_id"].iloc[0] == "ENST00000269305.8"


class TestValidation:
    def test_invalid_input_id(self, test_db):
        with GeneMapper(test_db) as gm:
            with pytest.raises(ValueError, match="Unknown id_type"):
                gm.map(["x"], input_id="bad_type", output_id="ncbi_gene_id")

    def test_invalid_output_id(self, test_db):
        with GeneMapper(test_db) as gm:
            with pytest.raises(ValueError, match="Unknown id_type"):
                gm.map(["x"], input_id="ensembl_gene_id", output_id="bad_type")

    def test_same_input_output(self, test_db):
        with GeneMapper(test_db) as gm:
            with pytest.raises(ValueError, match="must be different"):
                gm.map(["x"], input_id="ensembl_gene_id", output_id="ensembl_gene_id")

    def test_missing_db_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            GeneMapper(tmp_path / "nonexistent.db")


class TestContextManager:
    def test_context_manager(self, test_db):
        with GeneMapper(test_db) as gm:
            types = gm.list_id_types()
        assert len(types) > 0

    def test_empty_ids(self, test_db):
        with GeneMapper(test_db) as gm:
            mapped, unmapped = gm.map(
                [],
                input_id="ensembl_gene_id",
                output_id="ncbi_gene_id",
            )
        assert len(mapped) == 0
        assert len(unmapped) == 0
