"""
Build a geneXref database from an HGNC complete-set export.

Usage
-----
    python examples/build_database.py path/to/hgnc_complete_set.txt path/to/geneXref_database.tsv

The script demonstrates:
- Downloading or using a local copy of the HGNC complete set export
- Building a geneXref database with rebuild_database()
- Loading the new database and inspecting its contents
"""

import sys

from geneXref import geneXref


def main(hgnc_path: str, output_path: str) -> None:
    print(f"Building geneXref database from: {hgnc_path}")
    geneXref.rebuild_database(hgnc_path, output_path)
    print(f"Database written to: {output_path}")

    gx = geneXref(output_path)

    id_types = gx.list_id_types()
    print(f"\nAvailable identifier types ({len(id_types)}):")
    for id_type in id_types:
        print(f"  - {id_type}")

    print(f"\nTotal genes in database: {len(gx._db)}")
    print("\nSample entries:\n")
    print(gx._db.head(10).to_string(index=False))


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(
            f"Usage: python {sys.argv[0]} "
            f"path/to/hgnc_complete_set.txt path/to/geneXref_database.tsv"
        )
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
