# Changelog

All notable changes to this project will be documented in this file.

## 0.2.8 (unreleased)

- Argument validation throws `ArgumentError` or `DimensionMismatch`
  instead of a generic `ErrorException`. Code that catches
  `ErrorException` from the constructors needs updating.
- `Species`, `SNPSet` and trait definitions keep their multi-line display
  at the REPL (`MIME"text/plain"`) and gain a one-line `show`, so they
  print properly inside vectors and in string interpolation.
- `in(locus, ::LocusSet)` uses a binary search.
- `julia` compat is `1.10`, matching what CI tests.
- Tests run Aqua.jl; CI uploads coverage and no longer requests `x64` on
  macOS runners.

## 0.2.7 (2026-11-02)

- Documentation themes (JuliaBnG Sepia and Black). No change to the code.

## 0.2.6 (2026-09-29)

- Expanded docstrings for all exported types and functions. No change to
  the API.
- Documentation is deployed to <https://juliabng.github.io/BnGStructs/>.

## 0.2.5 (2026-08-14)

- Added variant maps, sorted locus sets, QTL architectures, and true breeding
  value calculations.

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
