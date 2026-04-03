import pytest
import pandas as pd


# ---------------------------------------------------------------------------
# Minimal pre-built database (as written by rebuild_database)
# ---------------------------------------------------------------------------

DB_ROWS = [
    # hgnc_id  gene_symbol  ncbi_gene_id  ensembl_gene_id   ucsc_id      refseq_accession  ensembl_transcript_id  uniprot_id
    ("1944",   "CHN2",      "1124",       "ENSG00000106069", "uc003szz.4", "NM_004067",     "ENST00000222792.11",  "P52757"),
    ("31456",  "C9orf153",  "389766",     "ENSG00000187753", "uc031teh.1", "NM_001010907",  "ENST00000339137.7",   "Q5TBE3"),
    ("17378",  "NAGPA",     "51172",      "ENSG00000103174", "uc002cyg.4", "NM_016256",     "ENST00000312251.8",   "Q9UK23"),
    ("6840",   "MAP2K1",    "5604",       "ENSG00000169032", "uc010bhq.4", "NM_001411065",  "ENST00000307102.10",  "Q02750"),
    # Row with several missing values (pseudogene)
    ("23122",  "SKA2P1",    "729012",     "ENSG00000232387", None,         "NG_027505",     None,                  None),
    # Row that shares a uniprot_id with another entry (many-to-many test)
    ("99999",  "DUPGENE",   "99999",      "ENSG00000999999", None,         None,            None,                  "P52757"),
]

COLUMNS = [
    "hgnc_id", "gene_symbol", "ncbi_gene_id", "ensembl_gene_id",
    "ucsc_id", "refseq_accession", "ensembl_transcript_id", "uniprot_id",
]


@pytest.fixture()
def db_path(tmp_path):
    """Write the minimal test database to a temp TSV and return its path."""
    df = pd.DataFrame(DB_ROWS, columns=COLUMNS)
    path = tmp_path / "geneXref_test.tsv"
    df.to_csv(path, sep="\t", index=False)
    return str(path)


# ---------------------------------------------------------------------------
# Minimal raw HGNC export (as downloaded from HGNC)
# ---------------------------------------------------------------------------

HGNC_ROWS = [
    # hgnc_id        symbol    entrez_id  ensembl_gene_id    ucsc_id       refseq_accession  uniprot_ids     mane_select
    ("HGNC:1944",   "CHN2",   "1124",    "ENSG00000106069", "uc003szz.4", "NM_004067",      "P52757",       "ENST00000222792.11|NM_004067.4"),
    ("HGNC:31456",  "C9orf153","389766", "ENSG00000187753", "uc031teh.1", "NM_001010907",   "Q5TBE3",       "ENST00000339137.7|NM_001276366.4"),
    ("HGNC:17378",  "NAGPA",  "51172",   "ENSG00000103174", "uc002cyg.4", "NM_016256",      "Q9UK23",       "ENST00000312251.8|NM_016256.4"),
    # Multiple uniprot_ids
    ("HGNC:6840",   "MAP2K1", "5604",    "ENSG00000169032", "uc010bhq.4", "NM_001411065",   "Q02750|A0A0A0MRZ7", "ENST00000307102.10|NM_002755.4"),
    # Missing optional fields
    ("HGNC:23122",  "SKA2P1", "729012",  "ENSG00000232387", None,         "NG_027505",      None,           None),
]

HGNC_COLUMNS = [
    "hgnc_id", "symbol", "entrez_id", "ensembl_gene_id",
    "ucsc_id", "refseq_accession", "uniprot_ids", "mane_select",
]


@pytest.fixture()
def hgnc_path(tmp_path):
    """Write minimal HGNC export to a temp TSV and return its path."""
    df = pd.DataFrame(HGNC_ROWS, columns=HGNC_COLUMNS)
    path = tmp_path / "hgnc_complete_set.txt"
    df.to_csv(path, sep="\t", index=False)
    return str(path)
