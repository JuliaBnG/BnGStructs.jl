"""
    struct SNPSet

Specification for a panel or designated set of SNP markers.

Useful for defining targeted SNP chips, low-density panels, or whole-genome sequence subsets
with minor allele frequency (MAF) criteria and exclusivity constraints.

# Fields
- `name::String`: Name or identifier of the SNP set (e.g. `"50kChip"`).
- `nlc::Int`: Number of loci in the set.
- `maf::Float64`: Minor allele frequency threshold (`0.0 <= maf < 0.5`).
- `exclusive::Bool`: Flag indicating whether loci in this set are exclusive (non-overlapping with other sets).

# Constructors
    SNPSet(name::String, nlc::Int; maf = 0.0, exclusive = false) -> SNPSet

# Examples
```julia
using BnGStructs

chip = SNPSet("50kChip", 50_000; maf = 0.05, exclusive = true)
chip.name      # "50kChip"
chip.nlc       # 50000
chip.maf       # 0.05
chip.exclusive # true
```
"""
struct SNPSet
    name::String
    nlc::Int
    maf::Float64
    exclusive::Bool
end

"""
    SNPSet(name::String, nlc::Int; maf = 0.0, exclusive = false) -> SNPSet

Construct a `SNPSet` with the given name, locus count, minor allele frequency (MAF) cutoff,
and exclusivity flag.

# Arguments
- `name::String`: Name or label for the SNP panel.
- `nlc::Int`: Number of loci (must be positive).

# Keywords
- `maf::Float64`: Minimum minor allele frequency threshold (default: `0.0`, must be in `[0.0, 0.5)`).
- `exclusive::Bool`: Whether the loci should be dedicated exclusively to this set (default: `false`).

# Returns
- `SNPSet`: Initialized SNP set specification.
"""
function SNPSet(name::String, nlc::Int; maf = 0.0, exclusive = false)
    0.0 ≤ maf < 0.5 || error("maf must be in [0.0, 0.5)")
    nlc > 0 || error("nlc must be positive")
    SNPSet(name, nlc, maf, exclusive)
end

"""
    Base.show(io::IO, snp::SNPSet)
Display the SNPSet object `snp` in a human-readable format.
"""
function Base.show(io::IO, snp::SNPSet)
    println(io, "                      SNPSet: ", snp.name)
    println(io, "              Number of loci: ", snp.nlc)
    println(io, "Minor allele frequency (MAF): ", snp.maf)
    println(io, "                   Exclusive: ", snp.exclusive)
end
