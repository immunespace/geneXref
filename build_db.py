#!/usr/bin/env python
"""CLI entry point for building the gene/transcript mapping database."""

import argparse

from geneXref.build import build_database


def main():
    parser = argparse.ArgumentParser(
        description="Build gene/transcript ID mapping SQLite database from Ensembl"
    )
    parser.add_argument(
        "--version", type=int, default=115,
        help="Ensembl release version (default: 115)",
    )
    parser.add_argument(
        "--output", type=str, default="geneXref.db",
        help="Output SQLite database path (default: geneXref.db)",
    )
    parser.add_argument(
        "--include-alt-haplotypes", action="store_true",
        help="Include genes on alternative haplotypes and patches (excluded by default)",
    )
    parser.add_argument(
        "--include-par-y", action="store_true",
        help="Include pseudoautosomal region (PAR) genes on chromosome Y (excluded by default)",
    )
    args = parser.parse_args()
    build_database(
        ensembl_version=args.version,
        output_path=args.output,
        filter_alt_haplotypes=not args.include_alt_haplotypes,
        filter_par_y=not args.include_par_y,
    )


if __name__ == "__main__":
    main()
