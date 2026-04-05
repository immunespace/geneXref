#!/usr/bin/env python
"""Example: build a geneXref SQLite database from the Ensembl public MySQL server.

Requires a network connection and the 'pymysql' package:
    pip install "geneXref[build]"

The resulting database is saved to ./geneXref.db by default.  Pass it to
GeneMapper(db_path=...) or move it to ~/.geneXref/geneXref.db so it is
picked up automatically.
"""

from geneXref.build import build_database

build_database(
    ensembl_version=115,        # Ensembl release to query
    output_path="geneXref.db",  # destination file
    filter_alt_haplotypes=True, # exclude alt-haplotype/patch sequences
    filter_par_y=True,          # exclude PAR_Y duplicates
)

print("Database built: geneXref.db")
