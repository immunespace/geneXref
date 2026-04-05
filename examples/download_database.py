#!/usr/bin/env python
"""Example: download a pre-built geneXref database from GitHub Releases.

The database is saved to ~/.geneXref/geneXref.db, where GeneMapper will
find it automatically.  No extra dependencies are required beyond the base
package install.
"""

from geneXref.download import download_db

db_path = download_db()
print(f"Ready to use: {db_path}")
