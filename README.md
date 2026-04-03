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

# initialize the geneXref object with the path to the database
geneXref = geneXref("path/to/geneXref_database.tsv")

idmap <- geneXref.map(["ENSG00000139618", "ENSG00000166710"],
                      input_id="ensembl_gene_id",
                      output_ids=["gene_symbol", "ncbi_gene_id"])
```

`idmap` will be a pandas DataFrame containing the mapping results, with columns for the input identifier and the requested output identifiers.

IDs that cannot be mapped will have `NaN` in the corresponding output columns, and a warning will be issued to indicate that no mapping was found.

IDs that map to multiple entries will also have `NaN` in the output columns, and a warning will be issued to indicate that multiple mappings were found.

Passing the argument `remove_unmapped=True` will remove rows with ambiguous or unmapped IDs from the output DataFrame.

### Listing available identifier types

```python
from geneXref import geneXref
geneXref = geneXref("path/to/geneXref_database.tsv")
geneXref.list_id_types()
```

### Rebuilding the database

1. Download the latest HGNC complete set export from the HGNC website and save it to your local machine.

```bash
wget https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt
```

2. Use the `rebuild_database` function to create a new geneXref database from the downloaded HGNC complete set file.
```python
from geneXref import geneXref
geneXref.rebuild_database("path/to/hgnc_complete_set.txt",
                          output_path="path/to/geneXref_database.tsv")
```
