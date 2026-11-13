# Changelog

All notable changes to this project will be documented in this file.

## 0.3.0 (unreleased)

Breaking release: the contracts that FisherWright, RelationshipMatrices,
Breeding and BnGsim rely on are now written down and enforced.

### Breaking

- `VariantMap.pos` is the genome-wide base-pair coordinate (chromosomes
  end to end), and the constructor requires it to be strictly increasing
  and `chr` to be non-decreasing. This was already how FisherWright and
  BnGsim used it; the documentation described a within-chromosome
  position and nothing was checked. Maps with within-chromosome positions
  must be built with `VariantMap(chr, chrpos, species)`. `vm[idxs]`
  requires increasing indices.
- `tbv` throws for a QTL architecture with non-zero dominance effects.
  It used to return the genotypic value under the name of a breeding
  value. Use `genotypic_value` or `breeding_value`.
- `VariantMap`, `TraitQTL` and `MultiTraitQTL` have a new last field
  (`morgan`, `intercept`), which changes their positional inner
  constructors. The documented constructors are unchanged.
- The unexported placeholder types `adTrait`, `adiTrait`, `dtTrait` and
  `ditTrait` are removed.
- `Species` is one concrete struct instead of an abstract type with ten
  identical subtypes. `Cattle(n)`, …, `GenericSpecies(name, n, lengths)`
  are constructor functions that return a `Species`. `x isa Species` and
  methods on `::Species` work as before; `x isa Cattle` and dispatch on
  `::Cattle` do not. `Species` has two new fields, `ends` and `morgan`.
- `name` and `M` are no longer exported, because they collide with
  ordinary variable names. Use `sp.name` and `sp.M`, or
  `BnGStructs.name(sp)` and `BnGStructs.M(sp)`.
- A genome beyond the `UInt32` coordinate range is refused when the
  species is constructed; `cbp` used to wrap around.

### New

- `chrpos`, `check_coordinates`, `VariantMap(pos, species)` and
  `VariantMap(chr, chrpos, species)`.
- Optional genetic positions (`morgan`, Morgans within the chromosome) in
  `VariantMap`.
- `genotypic_value` and `breeding_value` (average effects at stated
  reference frequencies).
- `intercept` in `TraitQTL` and `MultiTraitQTL`, so a base-population
  mean can be set without shifting the sampled effects.
- Genetic length per chromosome: `Species(...; morgan = ...)` and
  `map_length`. The default is `chromosome ./ M`, as before.
- Species without a population size: `Cattle()`, `GenericSpecies(name,
  lengths)`. The `nid` field and the constructors that take a size remain
  but are deprecated; `nid(sp)` throws when no size is set.
- `cbp(sp)` returns the stored chromosome ends instead of recomputing them.
- `select_individuals`, `select_loci`, `dosage` and `packed`.
- `write_haplotypes` and `read_haplotypes` for the `BNGH1` format.
- `trait_columns`: the record-column convention shared by BnGsim and
  Breeding.

### Not in this release

Removing `Species.nid` (once FisherWright and BnGsim take the population
size themselves), dropping `aTrait.QTL`, a single trait hierarchy, and sex
chromosomes are planned for later releases (see the JuliaBnG roadmap).

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
