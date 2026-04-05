# Gene and Transcript ID Mapping Database

A SQLite database for mapping gene and transcript identifiers across Ensembl, NCBI, HGNC, RefSeq, and UCSC. The [Ensembl public MySQL server](https://www.ensembl.org/info/data/mysql.html) is used as the source of truth.

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

- **ensg2all** — Joins `ensg` with `ensg2hgnc` and `ensg2ncbi` so you can look up a gene's HGNC and NCBI IDs in a single query.
- **enst2all** — Joins `enst` with `ensg`, `enst2refseq`, and `enst2ucsc` to show gene name, biotype, RefSeq, and UCSC IDs alongside each transcript.

## Usage

### Building the database

```bash
pip install -r requirements.txt
python build_db.py                    # Ensembl v115, writes ensg_enst_map.db
python build_db.py --version 116      # use a different Ensembl release
python build_db.py --output my.db     # custom output path
```

The script connects to `ensembldb.ensembl.org` (anonymous access, no credentials needed), fetches the relevant tables from `homo_sapiens_core_<version>_38`, removes multimappers, and writes the result to a local SQLite file.

### Querying

```sql
-- Gene info with all external IDs
SELECT * FROM ensg2all WHERE gene_name = 'BRCA1';

-- Transcript info with gene name and external IDs
SELECT * FROM enst2all WHERE gene_name = 'TP53';

-- Map an NCBI gene ID to Ensembl
SELECT * FROM ncbi2ensg WHERE ncbi_gene_id = '672';

-- Map a RefSeq transcript to Ensembl
SELECT * FROM refseq2enst WHERE refseq_mrna_id = 'NM_007294';
```

## How multimapper removal works

The raw Ensembl xref tables contain many-to-many relationships. To produce clean one-to-one mapping tables, the build script removes any identifier that maps to more than one partner. For example, when building `ensg2ncbi`, any NCBI gene ID that maps to multiple Ensembl gene IDs is dropped entirely. The forward and reverse tables are built independently, so they may contain slightly different subsets of the raw mappings.

## Source queries

The data is extracted from the Ensembl MySQL `homo_sapiens_core_<version>_38` schema using the external database IDs:

| Source | external_db_id |
|--------|----------------|
| NCBI | 1300 |
| HGNC | 1100 |
| RefSeq | 1801 |
| UCSC | 11000 |

See `build_db.py` for the full SQL queries.
