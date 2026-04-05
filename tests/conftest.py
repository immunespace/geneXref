"""Fixtures that create a small test database for GeneMapper tests."""

import sqlite3
import pytest

from geneXref.build import SQLITE_SCHEMA


@pytest.fixture()
def test_db(tmp_path):
    """Build a minimal SQLite database and return its path."""
    db_path = tmp_path / "test.db"
    conn = sqlite3.connect(str(db_path))
    conn.executescript(SQLITE_SCHEMA)

    # -- Core gene table --
    conn.executemany("INSERT INTO ensg VALUES (?, ?, ?, ?)", [
        ("ENSG00000141510", "protein_coding", "TP53", "tumor protein p53"),
        ("ENSG00000012048", "protein_coding", "BRCA1", "BRCA1 DNA repair"),
        ("ENSG00000146648", "protein_coding", "EGFR", "epidermal growth factor receptor"),
        ("ENSG00000099999", "lncRNA", "LNCRNA1", "test lncRNA"),
    ])

    # -- Transcript table --
    conn.executemany("INSERT INTO enst VALUES (?, ?, ?)", [
        ("ENST00000269305", "ENSG00000141510", "protein_coding"),
        ("ENST00000413465", "ENSG00000141510", "protein_coding"),
        ("ENST00000357654", "ENSG00000012048", "protein_coding"),
        ("ENST00000275493", "ENSG00000146648", "protein_coding"),
    ])

    # -- Gene mapping tables (1:1 after multimapper removal) --
    conn.executemany("INSERT INTO ensg2ncbi VALUES (?, ?)", [
        ("ENSG00000141510", "7157"),
        ("ENSG00000012048", "672"),
        ("ENSG00000146648", "1956"),
    ])
    conn.executemany("INSERT INTO ncbi2ensg VALUES (?, ?)", [
        ("7157", "ENSG00000141510"),
        ("672", "ENSG00000012048"),
        ("1956", "ENSG00000146648"),
    ])
    conn.executemany("INSERT INTO ensg2hgnc VALUES (?, ?)", [
        ("ENSG00000141510", "HGNC:11998"),
        ("ENSG00000012048", "HGNC:1100"),
        ("ENSG00000146648", "HGNC:3236"),
    ])
    conn.executemany("INSERT INTO hgnc2ensg VALUES (?, ?)", [
        ("HGNC:11998", "ENSG00000141510"),
        ("HGNC:1100", "ENSG00000012048"),
        ("HGNC:3236", "ENSG00000146648"),
    ])

    # -- Transcript mapping tables --
    conn.executemany("INSERT INTO enst2refseq VALUES (?, ?)", [
        ("ENST00000269305", "NM_000546"),
        ("ENST00000357654", "NM_007294"),
    ])
    conn.executemany("INSERT INTO refseq2enst VALUES (?, ?)", [
        ("NM_000546", "ENST00000269305"),
        ("NM_007294", "ENST00000357654"),
    ])
    conn.executemany("INSERT INTO enst2ucsc VALUES (?, ?)", [
        ("ENST00000269305", "uc002gig.2"),
        ("ENST00000275493", "uc003tpa.3"),
    ])
    conn.executemany("INSERT INTO ucsc2enst VALUES (?, ?)", [
        ("uc002gig.2", "ENST00000269305"),
        ("uc003tpa.3", "ENST00000275493"),
    ])

    conn.commit()
    conn.close()
    return db_path
