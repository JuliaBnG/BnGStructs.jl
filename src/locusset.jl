"""
    struct LocusSet
A named subset of variant loci identified by their 1-based indices into a `VariantMap`.

Commonly used to define:
- Visible SNP chip panels (e.g. 50k / 770k marker sets)
- Reference / background tracking loci
- Candidate causal regions
"""
struct LocusSet
    name::String
    loci::Vector{Int32}

    function LocusSet(name::AbstractString, loci::Vector{<:Integer})
        all(loci .> 0) || throw(ArgumentError("Locus indices must be positive integers"))
        issorted(loci) && allunique(loci) ||
            throw(ArgumentError("Locus indices must be sorted and unique"))
        new(String(name), Int32.(loci))
    end
end

Base.length(ls::LocusSet) = length(ls.loci)
Base.in(locus::Integer, ls::LocusSet) = locus in ls.loci

function Base.show(io::IO, ls::LocusSet)
    print(io, "LocusSet \"$(ls.name)\" with $(length(ls)) loci")
end
