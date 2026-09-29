"""
    struct LocusSet

A named subset of variant loci identified by their 1-based indices into a [`VariantMap`](@ref).

Commonly used to define:
- Visible SNP chip panels (e.g. 50k / 770k commercial marker arrays)
- Reference or background tracking loci
- Candidate causal or QTL candidate regions

# Fields
- `name::String`: Name or label of the locus set.
- `loci::Vector{Int32}`: Sorted, unique 1-based variant indices.

# Constructors
    LocusSet(name::AbstractString, loci::Vector{<:Integer}) -> LocusSet

Construct a `LocusSet`. Requires all indices in `loci` to be positive, strictly sorted in
ascending order, and unique.

# Methods
- `length(ls::LocusSet)`: Returns the number of loci in the set.
- `in(locus::Integer, ls::LocusSet)`: Membership test checking if a 1-based locus index is in the set.

# Examples
```julia
using BnGStructs

panel = LocusSet("50kChip", [1, 5, 12, 100])
length(panel) # 4
5 in panel    # true
10 in panel   # false
```
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
