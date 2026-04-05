"""
Build a SQLite database mapping gene/transcript IDs between Ensembl, NCBI, HGNC,
RefSeq, and UCSC using the Ensembl MySQL database as the source of truth.
"""

import argparse
import sqlite3
import pymysql


ENSEMBL_HOST = "ensembldb.ensembl.org"
ENSEMBL_PORT = 3306
ENSEMBL_USER = "anonymous"

# External DB IDs in Ensembl
EXTERNAL_DB = {
    "ncbi": 1300,
    "hgnc": 1100,
    "refseq": 1801,
    "ucsc": 11000,
}

# --------------------------------------------------------------------------- #
# Ensembl MySQL queries
# --------------------------------------------------------------------------- #

QUERY_ENSG = """
SELECT g.stable_id AS ensembl_gene_id,
       g.biotype AS gene_biotype,
       ga_name.`value` AS gene_name,
       g.description
  FROM gene g
  LEFT JOIN gene_attrib ga_name
    ON g.gene_id = ga_name.gene_id
   AND ga_name.attrib_type_id = 4
 WHERE g.stable_id NOT LIKE 'LRG%%'
"""

QUERY_ENST = """
SELECT g.stable_id AS ensembl_gene_id,
       t.stable_id AS ensembl_transcript_id,
       t.biotype AS transcript_biotype
  FROM gene g, transcript t
 WHERE g.gene_id = t.gene_id
   AND g.stable_id NOT LIKE 'LRG%%'
"""

QUERY_GENE_XREF = """
SELECT g.stable_id AS ensembl_gene_id,
       x.dbprimary_acc AS external_id
  FROM gene g, object_xref ox, xref x
 WHERE g.gene_id = ox.ensembl_id
   AND ox.xref_id = x.xref_id
   AND x.external_db_id = %s
   AND g.stable_id NOT LIKE 'LRG%%'
"""

QUERY_TRANSCRIPT_XREF = """
SELECT t.stable_id AS ensembl_transcript_id,
       x.dbprimary_acc AS external_id
  FROM transcript t, object_xref ox, xref x
 WHERE t.transcript_id = ox.ensembl_id
   AND ox.xref_id = x.xref_id
   AND x.external_db_id = %s
"""

# --------------------------------------------------------------------------- #
# SQLite schema
# --------------------------------------------------------------------------- #

SQLITE_SCHEMA = """
CREATE TABLE IF NOT EXISTS ensg (
    ensembl_gene_id TEXT PRIMARY KEY,
    gene_biotype TEXT,
    gene_name TEXT,
    description TEXT
);

CREATE TABLE IF NOT EXISTS enst (
    ensembl_transcript_id TEXT PRIMARY KEY,
    ensembl_gene_id TEXT NOT NULL,
    transcript_biotype TEXT,
    FOREIGN KEY (ensembl_gene_id) REFERENCES ensg(ensembl_gene_id)
);

CREATE TABLE IF NOT EXISTS ensg2ncbi (
    ensembl_gene_id TEXT PRIMARY KEY,
    ncbi_gene_id TEXT NOT NULL,
    FOREIGN KEY (ensembl_gene_id) REFERENCES ensg(ensembl_gene_id)
);

CREATE TABLE IF NOT EXISTS ensg2hgnc (
    ensembl_gene_id TEXT PRIMARY KEY,
    hgnc_id TEXT NOT NULL,
    FOREIGN KEY (ensembl_gene_id) REFERENCES ensg(ensembl_gene_id)
);

CREATE TABLE IF NOT EXISTS ncbi2ensg (
    ncbi_gene_id TEXT PRIMARY KEY,
    ensembl_gene_id TEXT NOT NULL,
    FOREIGN KEY (ensembl_gene_id) REFERENCES ensg(ensembl_gene_id)
);

CREATE TABLE IF NOT EXISTS hgnc2ensg (
    hgnc_id TEXT PRIMARY KEY,
    ensembl_gene_id TEXT NOT NULL,
    FOREIGN KEY (ensembl_gene_id) REFERENCES ensg(ensembl_gene_id)
);

CREATE TABLE IF NOT EXISTS enst2ucsc (
    ensembl_transcript_id TEXT PRIMARY KEY,
    ucsc_transcript_id TEXT NOT NULL,
    FOREIGN KEY (ensembl_transcript_id) REFERENCES enst(ensembl_transcript_id)
);

CREATE TABLE IF NOT EXISTS enst2refseq (
    ensembl_transcript_id TEXT PRIMARY KEY,
    refseq_mrna_id TEXT NOT NULL,
    FOREIGN KEY (ensembl_transcript_id) REFERENCES enst(ensembl_transcript_id)
);

CREATE TABLE IF NOT EXISTS ucsc2enst (
    ucsc_transcript_id TEXT PRIMARY KEY,
    ensembl_transcript_id TEXT NOT NULL,
    FOREIGN KEY (ensembl_transcript_id) REFERENCES enst(ensembl_transcript_id)
);

CREATE TABLE IF NOT EXISTS refseq2enst (
    refseq_mrna_id TEXT PRIMARY KEY,
    ensembl_transcript_id TEXT NOT NULL,
    FOREIGN KEY (ensembl_transcript_id) REFERENCES enst(ensembl_transcript_id)
);


CREATE VIEW IF NOT EXISTS ensg2all AS
SELECT g.*,
       h.hgnc_id,
       n.ncbi_gene_id
  FROM ensg g
  LEFT JOIN ensg2hgnc h ON g.ensembl_gene_id = h.ensembl_gene_id
  LEFT JOIN ensg2ncbi n ON g.ensembl_gene_id = n.ensembl_gene_id;

CREATE VIEW IF NOT EXISTS enst2all AS
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

def connect_ensembl(schema: str) -> pymysql.Connection:
    """Connect to the Ensembl public MySQL server."""
    return pymysql.connect(
        host=ENSEMBL_HOST,
        port=ENSEMBL_PORT,
        user=ENSEMBL_USER,
        database=schema,
        cursorclass=pymysql.cursors.Cursor,
    )


def fetch_all(conn: pymysql.Connection, query: str, params=None) -> list[tuple]:
    """Execute a query and return all rows."""
    with conn.cursor() as cur:
        cur.execute(query, params)
        return cur.fetchall()


def remove_multimappers(rows: list[tuple], key_col: int, val_col: int) -> list[tuple]:
    """
    Given rows of (key, value) pairs, remove any row where the value maps to
    more than one distinct key. Returns deduplicated rows.
    """
    # First, deduplicate identical pairs
    unique_pairs = list(set(rows))

    # Count how many distinct keys each value maps to
    val_to_keys: dict[str, set[str]] = {}
    for row in unique_pairs:
        val = row[val_col]
        key = row[key_col]
        val_to_keys.setdefault(val, set()).add(key)

    # Keep only values that map to exactly one key
    multi_vals = {v for v, keys in val_to_keys.items() if len(keys) > 1}

    return [row for row in unique_pairs if row[val_col] not in multi_vals]


def build_forward_map(rows: list[tuple]) -> list[tuple]:
    """
    Build the forward map (ensembl -> external): remove rows where the
    external ID maps to multiple ensembl IDs.
    key_col=0 is ensembl_id, val_col=1 is external_id.
    """
    return remove_multimappers(rows, key_col=0, val_col=1)


def build_reverse_map(rows: list[tuple]) -> list[tuple]:
    """
    Build the reverse map (external -> ensembl): remove rows where the
    ensembl ID maps to multiple external IDs.
    key_col=1 is external_id (becomes the key), val_col=0 is ensembl_id.
    """
    return remove_multimappers(rows, key_col=1, val_col=0)


# --------------------------------------------------------------------------- #
# Main build logic
# --------------------------------------------------------------------------- #

def build_database(ensembl_version: int, output_path: str):
    schema = f"homo_sapiens_core_{ensembl_version}_38"
    print(f"Connecting to {ENSEMBL_HOST}, schema: {schema}")
    mysql_conn = connect_ensembl(schema)

    # Fetch core tables
    print("Fetching gene table (ensg)...")
    ensg_rows = fetch_all(mysql_conn, QUERY_ENSG)
    print(f"  {len(ensg_rows)} genes")

    print("Fetching transcript table (enst)...")
    enst_rows = fetch_all(mysql_conn, QUERY_ENST)
    print(f"  {len(enst_rows)} transcripts")

    # Fetch gene xrefs
    print("Fetching NCBI gene xrefs...")
    ncbi_raw = fetch_all(mysql_conn, QUERY_GENE_XREF, (EXTERNAL_DB["ncbi"],))
    print(f"  {len(ncbi_raw)} raw mappings")

    print("Fetching HGNC gene xrefs...")
    hgnc_raw = fetch_all(mysql_conn, QUERY_GENE_XREF, (EXTERNAL_DB["hgnc"],))
    print(f"  {len(hgnc_raw)} raw mappings")

    # Fetch transcript xrefs
    print("Fetching UCSC transcript xrefs...")
    ucsc_raw = fetch_all(mysql_conn, QUERY_TRANSCRIPT_XREF, (EXTERNAL_DB["ucsc"],))
    print(f"  {len(ucsc_raw)} raw mappings")

    print("Fetching RefSeq transcript xrefs...")
    refseq_raw = fetch_all(mysql_conn, QUERY_TRANSCRIPT_XREF, (EXTERNAL_DB["refseq"],))
    print(f"  {len(refseq_raw)} raw mappings")

    mysql_conn.close()

    # Build clean mapping tables (remove multimappers)
    print("Removing multimappers...")

    ensg2ncbi = build_forward_map(ncbi_raw)
    ncbi2ensg = build_reverse_map(ncbi_raw)
    print(f"  ensg2ncbi: {len(ensg2ncbi)} (from {len(ncbi_raw)} raw)")
    print(f"  ncbi2ensg: {len(ncbi2ensg)} (from {len(ncbi_raw)} raw)")

    ensg2hgnc = build_forward_map(hgnc_raw)
    hgnc2ensg = build_reverse_map(hgnc_raw)
    print(f"  ensg2hgnc: {len(ensg2hgnc)} (from {len(hgnc_raw)} raw)")
    print(f"  hgnc2ensg: {len(hgnc2ensg)} (from {len(hgnc_raw)} raw)")

    enst2ucsc = build_forward_map(ucsc_raw)
    ucsc2enst = build_reverse_map(ucsc_raw)
    print(f"  enst2ucsc: {len(enst2ucsc)} (from {len(ucsc_raw)} raw)")
    print(f"  ucsc2enst: {len(ucsc2enst)} (from {len(ucsc_raw)} raw)")

    enst2refseq = build_forward_map(refseq_raw)
    refseq2enst = build_reverse_map(refseq_raw)
    print(f"  enst2refseq: {len(enst2refseq)} (from {len(refseq_raw)} raw)")
    print(f"  refseq2enst: {len(refseq2enst)} (from {len(refseq_raw)} raw)")

    # Write to SQLite
    print(f"Writing SQLite database to {output_path}...")
    sqlite_conn = sqlite3.connect(output_path)
    sqlite_conn.execute("PRAGMA foreign_keys = ON")
    sqlite_conn.executescript(SQLITE_SCHEMA)

    sqlite_conn.executemany(
        "INSERT OR IGNORE INTO ensg VALUES (?, ?, ?, ?)", ensg_rows
    )
    sqlite_conn.executemany(
        "INSERT OR IGNORE INTO enst VALUES (?, ?, ?)",
        [(t[1], t[0], t[2]) for t in enst_rows],  # reorder: enst, ensg, biotype
    )

    sqlite_conn.executemany(
        "INSERT OR IGNORE INTO ensg2ncbi VALUES (?, ?)", ensg2ncbi
    )
    sqlite_conn.executemany(
        "INSERT OR IGNORE INTO ensg2hgnc VALUES (?, ?)", ensg2hgnc
    )
    sqlite_conn.executemany(
        "INSERT OR IGNORE INTO ncbi2ensg VALUES (?, ?)",
        [(row[1], row[0]) for row in ncbi2ensg],
    )
    sqlite_conn.executemany(
        "INSERT OR IGNORE INTO hgnc2ensg VALUES (?, ?)",
        [(row[1], row[0]) for row in hgnc2ensg],
    )

    sqlite_conn.executemany(
        "INSERT OR IGNORE INTO enst2ucsc VALUES (?, ?)", enst2ucsc
    )
    sqlite_conn.executemany(
        "INSERT OR IGNORE INTO enst2refseq VALUES (?, ?)", enst2refseq
    )
    sqlite_conn.executemany(
        "INSERT OR IGNORE INTO ucsc2enst VALUES (?, ?)",
        [(row[1], row[0]) for row in ucsc2enst],
    )
    sqlite_conn.executemany(
        "INSERT OR IGNORE INTO refseq2enst VALUES (?, ?)",
        [(row[1], row[0]) for row in refseq2enst],
    )

    sqlite_conn.commit()

    # Print summary
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


def main():
    parser = argparse.ArgumentParser(
        description="Build gene/transcript ID mapping SQLite database from Ensembl"
    )
    parser.add_argument(
        "--version", type=int, default=115,
        help="Ensembl release version (default: 115)",
    )
    parser.add_argument(
        "--output", type=str, default="ensg_enst_map.db",
        help="Output SQLite database path (default: ensg_enst_map.db)",
    )
    args = parser.parse_args()
    build_database(args.version, args.output)


if __name__ == "__main__":
    main()
