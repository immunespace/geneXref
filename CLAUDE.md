# CLAUDE.md

## Project overview

geneXref is a lightweight Python package for mapping gene identifiers across
databases (Ensembl, NCBI, HGNC, UniProt, RefSeq, UCSC).  The source of truth
is the HGNC complete set export.  The package builds a flat TSV database from
that export and provides a `map()` method to translate between ID types.  A
pre-built database is bundled with the package so no manual download is needed.

## Repository layout

- `geneXref/geneXref.py` — single-module implementation (class `geneXref`)
- `tests/conftest.py` — shared fixtures (minimal test databases)
- `tests/test_geneXref.py` — pytest test suite
- `examples/` — runnable example scripts
- `geneXref/data/` — bundled pre-built database (shipped with the package)
- `data/` — HGNC source file (not committed; user-provided)

## Development

```bash
pip install -e ".[dev]"       # editable install with pytest
python -m pytest              # run the full test suite
```

All tests should pass before committing.  Tests create temporary databases via
fixtures in `conftest.py` — no external data files are needed to run them.

## Key conventions

- The database is a plain TSV loaded into a pandas DataFrame (`dtype=str`).
- All identifier values are stored as strings, never numeric types.
- The class name `geneXref` is intentionally lowercase (matches the package name).
- `map()` takes a single `output_id` (not a list).
- Ambiguous or missing mappings are set to `pd.NA`; a single summary
  `UserWarning` reports how many identifiers could not be mapped.
- Versioned Ensembl IDs are stripped before lookup but preserved in output.
- UniProt: only the first (primary) accession is kept per gene; see README
  "Design notes" for rationale.
- Database format is TSV (not pickle); see README "Design notes" for rationale.

## Testing patterns

- Test fixtures build minimal in-memory databases (5–6 rows) — see `conftest.py`.
- Tests are organized into classes by feature: `TestRebuildDatabase`,
  `TestListIdTypes`, `TestMap`.
- Warnings are tested with `pytest.warns(UserWarning, match=...)`.
- Validation errors are tested with `pytest.raises(ValueError, match=...)`.
