module BnGStructs

include("misc.jl")
# Structs for SNP genotypes, haplotype majored or ID majored.
include("haplotypes.jl")
include("genotypes.jl")
include("species.jl")
include("trait.jl")
# transformations between them
include("hap-gt.jl")
# SNP sets and locus sets
include("snpset.jl")
include("locusset.jl")
# Genomic variant map and QTL architectures
include("variantmap.jl")
include("qtl.jl")

export Haplotype, Genotype, hap2id, id2hap, Species
export Cat, Cattle, Chicken, Dog, GenericSpecies, Goat, Horse, Pig, Rabbit, Sheep
export AbstractTrait, Trait, aTrait, tTrait, SNPSet
export VariantMap, LocusSet, TraitQTL, MultiTraitQTL, tbv
export name, nid, chromosome, M, total_bp, cbp

end # module BnGStructs
