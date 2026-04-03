# geneXref

This package is intended to provide a simple interface to map gene identifiers
across different databases.  It relies on the mapping provided by the HGNC
complete set export, which is updated regularly and can be downloaded from the HGNC website.

## Installation

```bash
cd geneXref
pip install -e .
```

## Usage

### Mapping gene identifiers with the existing database

```python
from geneXref import geneXref

gx = geneXref("path/to/geneXref_database.tsv")

idmap = gx.map(["ENSG00000139618", "ENSG00000166710"],
               input_id="ensembl_gene_id",
               output_ids=["gene_symbol", "ncbi_gene_id"])
```

`idmap` will be a pandas DataFrame with one row per input identifier and columns
for the input identifier and each requested output identifier.

A mapping is set to `pd.NA` and a `UserWarning` is issued when:

- No row in the database matches the input identifier.
- More than one row matches the input identifier.
- The output value found is shared by a different row in the database
  (many-to-many relationship).

Passing `remove_unmapped=True` drops all rows that contain any `pd.NA` in the
output columns from the returned DataFrame.

Versioned Ensembl identifiers (e.g. `ENSG00000139618.14`, `ENST00000222792.11`)
are automatically mapped by stripping the version suffix before lookup.  The
original versioned value is preserved in the output DataFrame.

### Listing available identifier types

```python
from geneXref import geneXref

gx = geneXref("path/to/geneXref_database.tsv")
gx.list_id_types()
```

Returns a list of column names that can be used as `input_id` or `output_ids`
in `map()`.  The standard database produced by `rebuild_database` contains:

| Column | Description |
|---|---|
| `hgnc_id` | HGNC numeric identifier |
| `gene_symbol` | Approved HGNC gene symbol |
| `ncbi_gene_id` | NCBI (Entrez) gene identifier |
| `ensembl_gene_id` | Ensembl gene identifier |
| `ucsc_id` | UCSC genome browser identifier |
| `refseq_accession` | RefSeq accession |
| `ensembl_transcript_id` | Ensembl transcript identifier (MANE Select), version suffix stripped |
| `uniprot_id` | UniProt accession (primary entry) |

### Design notes

**Database format (TSV).**  The database is stored as a plain TSV file rather
than a binary format such as pickle.  While pickle loads roughly 4x faster, the
TSV loads in under 250 ms which is negligible for a one-time initialization cost.
TSV is human-inspectable, portable across pandas versions, and avoids the
security concerns associated with unpickling untrusted files.

**UniProt ID handling.**  The HGNC source data can contain multiple pipe-delimited
UniProt accessions per gene.  Only the first (primary) accession is retained.
In practice only ~0.3% of genes with a UniProt accession carry more than one ID,
so preserving all of them would add very little recall while introducing
many-to-many mapping complexity.

### Rebuilding the database

1. Download the latest HGNC complete set export from the HGNC website and save it to your local machine.

```bash
wget https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt
```

2. Use the `rebuild_database` static method to create a new geneXref database from the downloaded file.

```python
from geneXref import geneXref

geneXref.rebuild_database("path/to/hgnc_complete_set.txt",
                          output_path="path/to/geneXref_database.tsv")
```
