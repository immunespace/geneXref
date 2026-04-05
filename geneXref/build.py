"""Build the gene/transcript ID mapping SQLite database from Ensembl MySQL."""

from __future__ import annotations

import sqlite3

import pymysql


ENSEMBL_HOST = "ensembldb.ensembl.org"
ENSEMBL_PORT = 3306
ENSEMBL_USER = "anonymous"

EXTERNAL_DB = {
    "ncbi": 1300,
    "hgnc": 1100,
    "refseq": 1801,
    "ucsc": 11000,
}

# --------------------------------------------------------------------------- #
# Ensembl MySQL queries
# --------------------------------------------------------------------------- #

# Filter clause: exclude genes on alternative haplotypes and patches.
# These are identified by the 'non_ref' seq_region attribute
# (attrib_type_id = 16) and include MHC alt assemblies on chr6,
# novel/fix patches like HG1012_PATCH, etc.
_FILTER_ALT_HAPLOTYPES = """
   AND NOT EXISTS (
       SELECT 1 FROM seq_region_attrib sra
       WHERE sra.seq_region_id = sr.seq_region_id
       AND sra.attrib_type_id = 16
   )
"""

# Filter clause: exclude pseudoautosomal region (PAR) genes on chromosome Y.
# These are identical copies of their X-chromosome counterparts.
# The Y chromosome is not marked non_ref, so we filter by coordinate range:
#   PAR1: Y:10,001–2,781,479
#   PAR2: Y:56,887,903–57,217,415
_FILTER_PAR_Y = """
   AND NOT (sr.name = 'Y'
            AND ((g.seq_region_start >= 10001 AND g.seq_region_end <= 2781479)
              OR (g.seq_region_start >= 56887903 AND g.seq_region_end <= 57217415)))
"""


def _build_filter(filter_alt_haplotypes: bool, filter_par_y: bool) -> str:
    """Build the WHERE clause fragment for assembly filtering."""
    parts = []
    if filter_alt_haplotypes:
        parts.append(_FILTER_ALT_HAPLOTYPES)
    if filter_par_y:
        parts.append(_FILTER_PAR_Y)
    return "".join(parts)


def _query_ensg(assembly_filter: str) -> str:
    return f"""
SELECT g.stable_id AS ensembl_gene_id,
       g.biotype AS gene_biotype,
       ga_name.`value` AS gene_name,
       g.description
  FROM gene g
  JOIN seq_region sr ON g.seq_region_id = sr.seq_region_id
  LEFT JOIN gene_attrib ga_name
    ON g.gene_id = ga_name.gene_id
   AND ga_name.attrib_type_id = 4
 WHERE g.stable_id NOT LIKE 'LRG%%'
{assembly_filter}
"""


def _query_enst(assembly_filter: str) -> str:
    return f"""
SELECT g.stable_id AS ensembl_gene_id,
       t.stable_id AS ensembl_transcript_id,
       t.biotype AS transcript_biotype
  FROM gene g
  JOIN transcript t ON g.gene_id = t.gene_id
  JOIN seq_region sr ON g.seq_region_id = sr.seq_region_id
 WHERE g.stable_id NOT LIKE 'LRG%%'
{assembly_filter}
"""


def _query_gene_xref(assembly_filter: str) -> str:
    return f"""
SELECT g.stable_id AS ensembl_gene_id,
       x.dbprimary_acc AS external_id
  FROM gene g
  JOIN object_xref ox ON g.gene_id = ox.ensembl_id
  JOIN xref x ON ox.xref_id = x.xref_id
  JOIN seq_region sr ON g.seq_region_id = sr.seq_region_id
 WHERE x.external_db_id = %s
   AND g.stable_id NOT LIKE 'LRG%%'
{assembly_filter}
"""


def _query_transcript_xref(assembly_filter: str) -> str:
    return f"""
SELECT t.stable_id AS ensembl_transcript_id,
       x.dbprimary_acc AS external_id
  FROM transcript t
  JOIN gene g ON t.gene_id = g.gene_id
  JOIN object_xref ox ON t.transcript_id = ox.ensembl_id
  JOIN xref x ON ox.xref_id = x.xref_id
  JOIN seq_region sr ON g.seq_region_id = sr.seq_region_id
 WHERE x.external_db_id = %s
{assembly_filter}
"""


# --------------------------------------------------------------------------- #
# SQLite schema
# --------------------------------------------------------------------------- #

SQLITE_SCHEMA = """
DROP VIEW IF EXISTS enst2all;
DROP VIEW IF EXISTS ensg2all;
DROP TABLE IF EXISTS refseq2enst;
DROP TABLE IF EXISTS ucsc2enst;
DROP TABLE IF EXISTS enst2refseq;
DROP TABLE IF EXISTS enst2ucsc;
DROP TABLE IF EXISTS hgnc2ensg;
DROP TABLE IF EXISTS ncbi2ensg;
DROP TABLE IF EXISTS ensg2hgnc;
DROP TABLE IF EXISTS ensg2ncbi;
DROP TABLE IF EXISTS enst;
DROP TABLE IF EXISTS ensg;

CREATE TABLE ensg (
    ensembl_gene_id TEXT PRIMARY KEY,
    gene_biotype TEXT,
    gene_name TEXT,
    description TEXT
);

CREATE TABLE enst (
    ensembl_transcript_id TEXT PRIMARY KEY,
    ensembl_gene_id TEXT NOT NULL,
    transcript_biotype TEXT,
    FOREIGN KEY (ensembl_gene_id) REFERENCES ensg(ensembl_gene_id)
);

CREATE TABLE ensg2ncbi (
    ensembl_gene_id TEXT PRIMARY KEY,
    ncbi_gene_id TEXT NOT NULL,
    FOREIGN KEY (ensembl_gene_id) REFERENCES ensg(ensembl_gene_id)
);

CREATE TABLE ensg2hgnc (
    ensembl_gene_id TEXT PRIMARY KEY,
    hgnc_id TEXT NOT NULL,
    FOREIGN KEY (ensembl_gene_id) REFERENCES ensg(ensembl_gene_id)
);

CREATE TABLE ncbi2ensg (
    ncbi_gene_id TEXT PRIMARY KEY,
    ensembl_gene_id TEXT NOT NULL,
    FOREIGN KEY (ensembl_gene_id) REFERENCES ensg(ensembl_gene_id)
);

CREATE TABLE hgnc2ensg (
    hgnc_id TEXT PRIMARY KEY,
    ensembl_gene_id TEXT NOT NULL,
    FOREIGN KEY (ensembl_gene_id) REFERENCES ensg(ensembl_gene_id)
);

CREATE TABLE enst2ucsc (
    ensembl_transcript_id TEXT PRIMARY KEY,
    ucsc_transcript_id TEXT NOT NULL,
    FOREIGN KEY (ensembl_transcript_id) REFERENCES enst(ensembl_transcript_id)
);

CREATE TABLE enst2refseq (
    ensembl_transcript_id TEXT PRIMARY KEY,
    refseq_mrna_id TEXT NOT NULL,
    FOREIGN KEY (ensembl_transcript_id) REFERENCES enst(ensembl_transcript_id)
);

CREATE TABLE ucsc2enst (
    ucsc_transcript_id TEXT PRIMARY KEY,
    ensembl_transcript_id TEXT NOT NULL,
    FOREIGN KEY (ensembl_transcript_id) REFERENCES enst(ensembl_transcript_id)
);

CREATE TABLE refseq2enst (
    refseq_mrna_id TEXT PRIMARY KEY,
    ensembl_transcript_id TEXT NOT NULL,
    FOREIGN KEY (ensembl_transcript_id) REFERENCES enst(ensembl_transcript_id)
);

CREATE VIEW ensg2all AS
SELECT g.*,
       h.hgnc_id,
       n.ncbi_gene_id
  FROM ensg g
  LEFT JOIN ensg2hgnc h ON g.ensembl_gene_id = h.ensembl_gene_id
  LEFT JOIN ensg2ncbi n ON g.ensembl_gene_id = n.ensembl_gene_id;

CREATE VIEW enst2all AS
SELECT t.*,
       g.gene_name,
       g.gene_biotype,
       r.refseq_mrna_id,
       u.ucsc_transcript_id
  FROM enst t
  LEFT JOIN ensg g ON t.ensembl_gene_id = g.ensembl_gene_id
  LEFT JOIN enst2refseq r ON t.ensembl_transcript_id = r.ensembl_transcript_id
  LEFT JOIN enst2ucsc u ON t.ensembl_transcript_id = u.ensembl_transcript_id;
"""


# --------------------------------------------------------------------------- #
# Helper functions
# --------------------------------------------------------------------------- #

def _connect_ensembl(schema: str) -> pymysql.Connection:
    return pymysql.connect(
        host=ENSEMBL_HOST,
        port=ENSEMBL_PORT,
        user=ENSEMBL_USER,
        database=schema,
        cursorclass=pymysql.cursors.Cursor,
    )


def _fetch_all(conn: pymysql.Connection, query: str, params=None) -> list[tuple]:
    with conn.cursor() as cur:
        cur.execute(query, params)
        return cur.fetchall()


def _remove_multimappers(rows: list[tuple], key_col: int, val_col: int) -> list[tuple]:
    """Remove any row where the value maps to more than one distinct key."""
    unique_pairs = list(set(rows))
    val_to_keys: dict[str, set[str]] = {}
    for row in unique_pairs:
        val_to_keys.setdefault(row[val_col], set()).add(row[key_col])
    multi_vals = {v for v, keys in val_to_keys.items() if len(keys) > 1}
    return [row for row in unique_pairs if row[val_col] not in multi_vals]


def _build_forward_map(rows: list[tuple]) -> list[tuple]:
    return _remove_multimappers(rows, key_col=0, val_col=1)


def _build_reverse_map(rows: list[tuple]) -> list[tuple]:
    return _remove_multimappers(rows, key_col=1, val_col=0)


# --------------------------------------------------------------------------- #
# Public build function
# --------------------------------------------------------------------------- #

def build_database(
    ensembl_version: int = 115,
    output_path: str = "geneXref.db",
    filter_alt_haplotypes: bool = True,
    filter_par_y: bool = True,
):
    """Fetch data from Ensembl MySQL and write a local SQLite database.

    Parameters
    ----------
    ensembl_version : int
        Ensembl release number (e.g. 115).
    output_path : str
        Destination path for the SQLite file.
    filter_alt_haplotypes : bool
        If True (default), exclude genes on alternative haplotypes and patches.
    filter_par_y : bool
        If True (default), exclude PAR genes on chromosome Y.
    """
    assembly_filter = _build_filter(filter_alt_haplotypes, filter_par_y)

    schema = f"homo_sapiens_core_{ensembl_version}_38"
    print(f"Connecting to {ENSEMBL_HOST}, schema: {schema}")
    mysql_conn = _connect_ensembl(schema)

    filters = []
    if filter_alt_haplotypes:
        filters.append("alt haplotypes/patches")
    if filter_par_y:
        filters.append("PAR_Y")
    if filters:
        print(f"Filtering: {', '.join(filters)}")
    else:
        print("No assembly filtering applied")

    query_ensg = _query_ensg(assembly_filter)
    query_enst = _query_enst(assembly_filter)
    query_gene_xref = _query_gene_xref(assembly_filter)
    query_transcript_xref = _query_transcript_xref(assembly_filter)

    print("Fetching gene table (ensg)...")
    ensg_rows = _fetch_all(mysql_conn, query_ensg)
    print(f"  {len(ensg_rows)} genes")

    print("Fetching transcript table (enst)...")
    enst_rows = _fetch_all(mysql_conn, query_enst)
    print(f"  {len(enst_rows)} transcripts")

    print("Fetching NCBI gene xrefs...")
    ncbi_raw = _fetch_all(mysql_conn, query_gene_xref, (EXTERNAL_DB["ncbi"],))
    print(f"  {len(ncbi_raw)} raw mappings")

    print("Fetching HGNC gene xrefs...")
    hgnc_raw = _fetch_all(mysql_conn, query_gene_xref, (EXTERNAL_DB["hgnc"],))
    print(f"  {len(hgnc_raw)} raw mappings")

    print("Fetching UCSC transcript xrefs...")
    ucsc_raw = _fetch_all(mysql_conn, query_transcript_xref, (EXTERNAL_DB["ucsc"],))
    print(f"  {len(ucsc_raw)} raw mappings")

    print("Fetching RefSeq transcript xrefs...")
    refseq_raw = _fetch_all(mysql_conn, query_transcript_xref, (EXTERNAL_DB["refseq"],))
    print(f"  {len(refseq_raw)} raw mappings")

    mysql_conn.close()

    print("Removing multimappers...")
    ensg2ncbi = _build_forward_map(ncbi_raw)
    ncbi2ensg = _build_reverse_map(ncbi_raw)
    print(f"  ensg2ncbi: {len(ensg2ncbi)} (from {len(ncbi_raw)} raw)")
    print(f"  ncbi2ensg: {len(ncbi2ensg)} (from {len(ncbi_raw)} raw)")

    ensg2hgnc = _build_forward_map(hgnc_raw)
    hgnc2ensg = _build_reverse_map(hgnc_raw)
    print(f"  ensg2hgnc: {len(ensg2hgnc)} (from {len(hgnc_raw)} raw)")
    print(f"  hgnc2ensg: {len(hgnc2ensg)} (from {len(hgnc_raw)} raw)")

    enst2ucsc = _build_forward_map(ucsc_raw)
    ucsc2enst = _build_reverse_map(ucsc_raw)
    print(f"  enst2ucsc: {len(enst2ucsc)} (from {len(ucsc_raw)} raw)")
    print(f"  ucsc2enst: {len(ucsc2enst)} (from {len(ucsc_raw)} raw)")

    enst2refseq = _build_forward_map(refseq_raw)
    refseq2enst = _build_reverse_map(refseq_raw)
    print(f"  enst2refseq: {len(enst2refseq)} (from {len(refseq_raw)} raw)")
    print(f"  refseq2enst: {len(refseq2enst)} (from {len(refseq_raw)} raw)")

    print(f"Writing SQLite database to {output_path}...")
    sqlite_conn = sqlite3.connect(output_path)
    sqlite_conn.execute("PRAGMA foreign_keys = ON")
    sqlite_conn.executescript(SQLITE_SCHEMA)

    sqlite_conn.executemany(
        "INSERT OR IGNORE INTO ensg VALUES (?, ?, ?, ?)", ensg_rows
    )
    sqlite_conn.executemany(
        "INSERT OR IGNORE INTO enst VALUES (?, ?, ?)",
        [(t[1], t[0], t[2]) for t in enst_rows],
    )
    sqlite_conn.executemany("INSERT OR IGNORE INTO ensg2ncbi VALUES (?, ?)", ensg2ncbi)
    sqlite_conn.executemany("INSERT OR IGNORE INTO ensg2hgnc VALUES (?, ?)", ensg2hgnc)
    sqlite_conn.executemany(
        "INSERT OR IGNORE INTO ncbi2ensg VALUES (?, ?)",
        [(r[1], r[0]) for r in ncbi2ensg],
    )
    sqlite_conn.executemany(
        "INSERT OR IGNORE INTO hgnc2ensg VALUES (?, ?)",
        [(r[1], r[0]) for r in hgnc2ensg],
    )
    sqlite_conn.executemany("INSERT OR IGNORE INTO enst2ucsc VALUES (?, ?)", enst2ucsc)
    sqlite_conn.executemany("INSERT OR IGNORE INTO enst2refseq VALUES (?, ?)", enst2refseq)
    sqlite_conn.executemany(
        "INSERT OR IGNORE INTO ucsc2enst VALUES (?, ?)",
        [(r[1], r[0]) for r in ucsc2enst],
    )
    sqlite_conn.executemany(
        "INSERT OR IGNORE INTO refseq2enst VALUES (?, ?)",
        [(r[1], r[0]) for r in refseq2enst],
    )

    sqlite_conn.commit()

    print("\nDatabase summary:")
    for table in [
        "ensg", "enst",
        "ensg2ncbi", "ensg2hgnc", "ncbi2ensg", "hgnc2ensg",
        "enst2ucsc", "enst2refseq", "ucsc2enst", "refseq2enst",
    ]:
        count = sqlite_conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
        print(f"  {table}: {count} rows")

    sqlite_conn.close()
    print("Done.")
