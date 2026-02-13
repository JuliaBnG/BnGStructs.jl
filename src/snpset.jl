"""
    struct SNPSet
        name::String
        nlc::Int
        maf::Float64
        exclusive::Bool
    end
"""
struct SNPSet
    name::String
    nlc::Int
    maf::Float64
    exclusive::Bool
end

"""
    SNPSet(name::String, nlc::Int; maf = 0.0, exclusive = false)
Create a SNPSet object with the given name, number of loci, minor allele
frequency (MAF), and exclusivity flag.
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
