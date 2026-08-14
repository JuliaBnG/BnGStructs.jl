# Changelog

All notable changes to this project will be documented in this file.

## 0.2.4 (2026-08-14)

- Added package documentation with a usage guide and API reference.
- Made species and trait storage type-stable and added validated species accessors.
- Added `AbstractMatrix` support and concise display output for haplotypes and genotypes.
- Hardened threshold-trait input validation and preserved caller-provided weights.

## 0.2.2 (2026-02-13)

- Bumped package version to 0.2.2
- Added `snpset.jl` to define SNP sets.
- Enabled `trait.jl` for trait definitions.

## 0.2.0 (2025-11-12)

- Bumped package version to 0.2.0.
- Added `changelog.md` with release notes.
- Added `Species.jl` for various normal species. Sex and mitochondrial
  chromosomes are ignored.
- Tests: existing test suite passes (see `test/runtests.jl`).

Notes:
- No breaking API changes expected compared to 0.1.3 (assumption based on
  current sources and tests).
