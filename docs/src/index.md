# BnGStructs.jl

BnGStructs.jl provides compact data structures for SNP haplotypes, genotypes,
species metadata, trait definitions, and SNP sets used in breeding and
population-genetics workflows.

## Installation

Install BnGStructs from the Julia General registry:

```julia
using Pkg
Pkg.add("BnGStructs")
```

## Quick start

```julia
using BnGStructs

# 1,000 loci and 100 diploid individuals (200 haplotypes).
haplotypes = Haplotype(1_000, 200)
genotypes = hap2id(haplotypes)

# Built-in species metadata and trait definitions.
cattle = Cattle(100)
milk_yield = Trait("MilkYield"; h²=0.3, σₐ=2.5)
disease_risk = Trait("DiseaseRisk", [0.7, 0.2, 0.1]; h²=0.2)
```

`Haplotype` stores loci in rows and haplotypes in columns. `Genotype` stores
individuals in rows and alleles in columns. Both are `AbstractMatrix{Bool}`
implementations; `size`, scalar indexing, and scalar assignment operate on the
logical matrix without exposing storage padding.

## Species metadata

`Cattle`, `Sheep`, `Chicken`, `Pig`, `Goat`, `Horse`, `Rabbit`, `Cat`, and
`Dog` supply curated autosomal chromosome lengths. Use `GenericSpecies` for a
custom genome:

```julia
species = GenericSpecies("Example", 500, [10_000_000, 20_000_000])
total_bp(species)
cbp(species)
```

The exported `name`, `nid`, `chromosome`, and `M` accessors provide the core
metadata. `total_bp` returns total genome length, while `cbp` returns cumulative
chromosome endpoints.

## Traits

`Trait(name; ...)` creates a continuous additive trait. Threshold traits use
relative category weights:

```julia
trait = Trait("DiseaseStatus", [0.95, 0.05]; h²=0.15)
```

Weights must be finite, positive, and contain at least two categories. They are
normalized internally without modifying the input vector.

## Variant maps and QTL architectures

`VariantMap` stores per-locus metadata, while `LocusSet` names a sorted, unique
subset of its one-based locus indices. `TraitQTL` and `MultiTraitQTL` represent
additive and optional dominance effects; `tbv` calculates true breeding values
from a `Haplotype` or `Genotype`:

```julia
map = VariantMap([1, 1], [100, 200], ['A', 'C'], ['G', 'T'])
panel = LocusSet("Example panel", [1, 2])
qtl = TraitQTL("Example trait", [1], [0.5])
values = tbv(haplotypes, qtl)
```
