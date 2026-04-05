#!/usr/bin/env python
"""Example: convert Ensembl gene IDs to HGNC gene names.

Requires the geneXref database to be present (run download_database.py first).
"""

from geneXref import GeneMapper

gene_ids = [
    "ENSG00000141510",  # TP53
    "ENSG00000012048",  # BRCA1
    "ENSG00000139618",  # BRCA2
    "ENSG00000157764",  # BRAF
    "ENSG00000999999",  # fictitious — will appear in unmapped
]

with GeneMapper() as mapper:
    mapped, unmapped = mapper.map(
        ids=gene_ids,
        input_id="ensembl_gene_id",
        output_id="gene_name",
    )

print("Mapped:")
print(mapped.to_string(index=False))

if not unmapped.empty:
    print("\nUnmapped:")
    print(unmapped.to_string(index=False))
