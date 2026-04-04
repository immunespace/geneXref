"""
Map a list of gene symbols to Ensembl gene IDs.

Usage
-----
    python examples/map_gene_symbols.py

The script demonstrates:
- Basic identifier mapping with geneXref.map()
- Inspecting mapped results and unmapped-ID reasons
"""

import sys
import warnings

from geneXref import geneXref


def main() -> None:
    gx = geneXref()

    gene_symbols = [
        "TP53",
        "BRCA1",
        "EGFR",
        "THIS_GENE_DOES_NOT_EXIST",
    ]

    print("=== Mapping gene symbols ===\n")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        mapped, unmapped = gx.map(
            gene_symbols,
            input_id="gene_symbol",
            output_id="ensembl_gene_id",
        )

    for w in caught:
        print(f"Warning: {w.message}")

    print(f"\nMapped ({len(mapped)}):\n")
    print(mapped.to_string(index=False))

    if len(unmapped) > 0:
        print(f"\nUnmapped ({len(unmapped)}):\n")
        print(unmapped["reason"].value_counts().to_string())


if __name__ == "__main__":
    if len(sys.argv) != 1:
        print(f"Usage: python {sys.argv[0]}")
        sys.exit(1)
    main()
