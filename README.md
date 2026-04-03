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
| `ensembl_transcript_id` | Ensembl transcript identifier (MANE Select) |
| `uniprot_id` | UniProt accession (primary entry) |

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
