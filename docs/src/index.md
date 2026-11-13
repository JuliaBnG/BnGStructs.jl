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
subset of its one-based locus indices.

Positions in a `VariantMap` are **genome-wide** base-pair coordinates:
chromosomes are laid end to end in the order of the species, and `pos` is
strictly increasing. This is the coordinate FisherWright simulates on and
BnGsim's meiosis searches. Build a map from within-chromosome positions
with a species, and recover them with `chrpos`:

```julia
species = GenericSpecies("Example", 500, [10_000, 20_000])
map = VariantMap([1, 1, 2], [100, 200, 150], species)   # chr, position on the chromosome
map.pos                                                # UInt32[100, 200, 10150]
chrpos(map, species)                                   # UInt32[100, 200, 150]
check_coordinates(map, species)                        # true
panel = LocusSet("Example panel", [1, 2])
```

`TraitQTL` and `MultiTraitQTL` hold additive effects, optional dominance
effects and an intercept. Three functions turn genotypes into values, and
they differ once there is dominance:

| Function | Value | Use |
| :--- | :--- | :--- |
| `genotypic_value` | ``G = \mu + \sum z a + \sum_{\text{het}} d`` | what a phenotype is built on |
| `tbv` | ``\mu + \sum z a`` | additive architectures only; throws if there is dominance |
| `breeding_value` | ``A = \sum (z - 2p)\alpha``, ``\alpha = a + d(1 - 2p)`` | what an EBV estimates; depends on the reference frequencies `p` |

```julia
qtl = TraitQTL("Example trait", [1], [0.5]; intercept = 10.0)
values = tbv(haplotypes, qtl)
```

## Subsets, dosages and files

`select_individuals` and `select_loci` return new `Haplotype`s, `dosage`
returns allele counts, and `packed` gives other packages the bit storage.
`write_haplotypes` and `read_haplotypes` use the `BNGH1` format, which is
also what the msprime founder scripts write. `trait_columns` defines the
column names (`tbv_<name>`, `ft_<name>`, `ebv_<name>`) under which BnGsim
and Breeding exchange records.
