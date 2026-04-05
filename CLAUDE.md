# CLAUDE.md

## Project overview

geneXref maps gene and transcript identifiers across Ensembl, NCBI, HGNC, RefSeq, and UCSC. It uses a SQLite database built from the Ensembl public MySQL server as the source of truth.

## Package structure

- `geneXref/` — the installable Python package
  - `mapper.py` — `GeneMapper` class: the main API (`map()`, `list_id_types()`)
  - `build.py` — database builder: queries Ensembl MySQL, removes multimappers, writes SQLite
  - `download.py` — downloads pre-built database from GitHub Releases to `~/.geneXref/`
  - `__init__.py` — lazy imports to avoid pulling in pandas when only building
- `build_db.py` — CLI entry point for building the database
- `tests/` — pytest suite with an in-memory SQLite fixture (no network needed)

## Key commands

```bash
pip install -e ".[dev]"       # install with dev dependencies
python -m pytest tests/ -v    # run tests
python build_db.py             # build database from Ensembl (requires network)
python -m geneXref.download   # download pre-built database from GitHub Releases
```

## Architecture notes

- The `__init__.py` uses `__getattr__` for lazy imports so that `from geneXref.build import build_database` doesn't pull in pandas.
- Build queries join `seq_region` and filter out non-primary-assembly genes (alt haplotypes via `non_ref` attrib, PAR_Y by coordinate range). These filters are controlled by `--include-alt-haplotypes` and `--include-par-y` flags.
- The mapper resolves queries through SQLite views (`ensg2all`, `enst2all`). Cross-level mappings (e.g. gene_name → refseq_mrna_id) join both views automatically.
- Default DB path resolution: `~/.geneXref/geneXref.db` → bundled `data/geneXref.db` → error with instructions.

## Testing

Tests use a small fixture database created in `tests/conftest.py` from the `SQLITE_SCHEMA` constant — no Ensembl connection needed. The fixture covers gene-level, transcript-level, and cross-level mappings plus edge cases (not_found, duplicates, version stripping, validation).
