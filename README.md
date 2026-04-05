# geneXref

A Python package for mapping gene and transcript identifiers across Ensembl, NCBI, HGNC, RefSeq, and UCSC. Backed by a SQLite database built from the [Ensembl public MySQL server](https://www.ensembl.org/info/data/mysql.html).

## Installation

```bash
pip install -e .            # core package (mapper only)
pip install -e ".[build]"   # includes pymysql for building the database
pip install -e ".[dev]"     # includes pymysql + pytest
```

## Quick start

### Downloading a pre-built database

The easiest way to get started is to download a pre-built database from [GitHub Releases](https://github.com/immunespace/geneXref/releases):

```bash
python -m geneXref.download
```

This saves `geneXref.db` to `~/.geneXref/`, where `GeneMapper()` will find it automatically.

### Building the database from Ensembl

```bash
python build_db.py                          # Ensembl v115 -> geneXref.db
python build_db.py --version 116            # use a different Ensembl release
python build_db.py --output my.db           # custom output path
python build_db.py --include-alt-haplotypes # include alt haplotypes and patches
python build_db.py --include-par-y          # include PAR genes on chromosome Y
```

The script connects to `ensembldb.ensembl.org` (anonymous, no credentials needed), fetches gene and transcript cross-references from `homo_sapiens_core_<version>_38`, removes multimappers, and writes the result to a local SQLite file.

By default, genes on alternative haplotypes, patches, and PAR regions of chromosome Y are excluded (see [Primary assembly filtering](#primary-assembly-filtering)). Use `--include-alt-haplotypes` and/or `--include-par-y` to include them.

### Using the Python API

```python
from geneXref import GeneMapper

# Uses ~/.geneXref/geneXref.db (downloaded) or bundled DB automatically
gm = GeneMapper()

# Or specify a path explicitly
gm = GeneMapper("geneXref.db")

# List available identifier types
gm.list_id_types()
# ['ensembl_gene_id', 'ensembl_transcript_id', 'gene_name',
#  'hgnc_id', 'ncbi_gene_id', 'refseq_mrna_id', 'ucsc_transcript_id']

# Map gene symbols to Ensembl gene IDs
mapped, unmapped = gm.map(
    ["TP53", "BRCA1", "EGFR"],
    input_id="gene_name",
    output_id="ensembl_gene_id",
)

# Map NCBI gene IDs to HGNC (chains through Ensembl automatically)
mapped, unmapped = gm.map(
    ["7157", "672"],
    input_id="ncbi_gene_id",
    output_id="hgnc_id",
)

# Cross-level mapping: gene name to RefSeq transcript
mapped, unmapped = gm.map(
    ["TP53"],
    input_id="gene_name",
    output_id="refseq_mrna_id",
)

# Ensembl version suffixes are stripped automatically
mapped, unmapped = gm.map(
    ["ENSG00000141510.14"],
    input_id="ensembl_gene_id",
    output_id="ncbi_gene_id",
)
```

`map()` returns a tuple of two DataFrames:

- **mapped** — successfully converted identifiers with columns `[input_id, output_id]`
- **unmapped** — failed conversions with columns `["identifier", "reason"]` where reason is one of `"not_found"`, `"duplicate_input"`, or `"duplicate_output"`

A `UserWarning` is issued when any identifiers cannot be mapped.

### Querying the database directly

```sql
SELECT * FROM ensg2all WHERE gene_name = 'BRCA1';
SELECT * FROM enst2all WHERE gene_name = 'TP53';
SELECT * FROM ncbi2ensg WHERE ncbi_gene_id = '672';
SELECT * FROM refseq2enst WHERE refseq_mrna_id = 'NM_007294';
```

## Schema

### Core tables

- **ensg** — One row per Ensembl gene. Columns: `ensembl_gene_id` (PK), `gene_biotype`, `gene_name`, `description`.
- **enst** — One row per Ensembl transcript. Columns: `ensembl_transcript_id` (PK), `ensembl_gene_id` (FK to `ensg`), `transcript_biotype`.

### Gene mapping tables

Each table contains only one-to-one mappings; any identifier that maps to multiple partners is removed.

| Table | Primary key | Maps to |
|-------|-------------|---------|
| ensg2ncbi | ensembl_gene_id | ncbi_gene_id |
| ensg2hgnc | ensembl_gene_id | hgnc_id |
| ncbi2ensg | ncbi_gene_id | ensembl_gene_id |
| hgnc2ensg | hgnc_id | ensembl_gene_id |

### Transcript mapping tables

Same one-to-one constraint as the gene mapping tables.

| Table | Primary key | Maps to |
|-------|-------------|---------|
| enst2ucsc | ensembl_transcript_id | ucsc_transcript_id |
| enst2refseq | ensembl_transcript_id | refseq_mrna_id |
| ucsc2enst | ucsc_transcript_id | ensembl_transcript_id |
| refseq2enst | refseq_mrna_id | ensembl_transcript_id |

### Convenience views

- **ensg2all** — Joins `ensg` with `ensg2hgnc` and `ensg2ncbi` for a unified gene-level view.
- **enst2all** — Joins `enst` with `ensg`, `enst2refseq`, and `enst2ucsc` for a unified transcript-level view.

## Primary assembly filtering

The GRCh38 assembly includes alternative haplotypes (e.g. MHC region on chromosome 6), fix/novel patches, and pseudoautosomal regions (PAR) on the Y chromosome. Ensembl creates separate gene models for each of these, but they are redundant copies of primary-assembly genes sharing the same gene symbol, HGNC ID, and NCBI gene ID. If included, they produce spurious `duplicate_output` / `duplicate_input` results.

By default, the build queries restrict to the primary assembly by excluding:

1. **Alternative haplotypes and patches** — identified by the `non_ref` attribute (`seq_region_attrib.attrib_type_id = 16`) in the Ensembl schema. This covers alt assemblies like `HSCHR6_MHC_COX_CTG1` and patches like `HG1012_PATCH`. Use `--include-alt-haplotypes` to include these.

2. **PAR genes on chromosome Y** — the Y chromosome is not marked `non_ref`, so PAR regions are filtered by coordinate range (PAR1: Y:10,001–2,781,479; PAR2: Y:56,887,903–57,217,415). The X-chromosome copies are retained. Use `--include-par-y` to include these.

This matches the convention used by HGNC, which tracks only one entry per gene on the primary assembly.

## How multimapper removal works

The raw Ensembl xref tables contain many-to-many relationships. To produce clean one-to-one mapping tables, the build script removes any identifier that maps to more than one partner. For example, when building `ensg2ncbi`, any NCBI gene ID that maps to multiple Ensembl gene IDs is dropped entirely. The forward and reverse tables are built independently, so they may contain slightly different subsets of the raw mappings.

## Source queries

Data is extracted from the Ensembl MySQL `homo_sapiens_core_<version>_38` schema using these external database IDs:

| Source | external_db_id |
|--------|----------------|
| NCBI | 1300 |
| HGNC | 1100 |
| RefSeq | 1801 |
| UCSC | 11000 |

See [geneXref/build.py](geneXref/build.py) for the full SQL queries.
