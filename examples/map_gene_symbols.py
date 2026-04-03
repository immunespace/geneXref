"""
Map a list of gene symbols to Ensembl gene IDs and NCBI gene IDs.

Usage
-----
    python examples/map_gene_symbols.py path/to/geneXref_database.tsv

The script demonstrates:
- Basic identifier mapping with geneXref.map()
- Handling of unmapped / ambiguous identifiers via warnings
- Filtering to confidently mapped rows with remove_unmapped=True
"""

import sys
import warnings

from geneXref import geneXref


def main(db_path: str) -> None:
    gx = geneXref(db_path)

    gene_symbols = [
        "TP53",
        "BRCA1",
        "EGFR",
        "THIS_GENE_DOES_NOT_EXIST",
    ]

    print("=== Mapping gene symbols (all results, including unmapped) ===\n")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = gx.map(
            gene_symbols,
            input_id="gene_symbol",
            output_ids=["ensembl_gene_id", "ncbi_gene_id"],
        )

    for w in caught:
        print(f"Warning: {w.message}")

    print()
    print(result.to_string(index=False))

    print("\n=== Mapped results only (remove_unmapped=True) ===\n")
    with warnings.catch_warnings(record=True):
        warnings.simplefilter("always")
        result_clean = gx.map(
            gene_symbols,
            input_id="gene_symbol",
            output_ids=["ensembl_gene_id", "ncbi_gene_id"],
            remove_unmapped=True,
        )

    print(result_clean.to_string(index=False))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: python {sys.argv[0]} path/to/geneXref_database.tsv")
        sys.exit(1)
    main(sys.argv[1])
