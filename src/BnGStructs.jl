module BnGStructs

include("misc.jl")
# Structs for SNP genotypes, haplotype majored or ID majored.
include("haplotypes.jl")
include("genotypes.jl")
include("species.jl")
include("trait.jl")
# transformations between them
include("hap-gt.jl")

export Haplotype, Genotype, hap2id, id2hap, Species, Trait

end # module BnGStructs
