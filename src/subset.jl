"""
    packed(hap::Haplotype) -> BitMatrix
    packed(gt::Genotype) -> BitMatrix

The bit storage of `hap` or `gt`, not a copy. Its first dimension is padded
to a multiple of 64 (so that each column starts on a 64-bit word and
columns can be written from different threads), and the padding bits are
zero. Rows `1:nlc` (or `1:nid`) hold the data. This is the supported way
for other packages to reach the words; treat the padding as read-only.
"""
packed(hap::Haplotype) = hap.gt
packed(gt::Genotype) = gt.gt

"""
    select_individuals(hap::Haplotype, ids) -> Haplotype

A new `Haplotype` with the two haplotypes of each individual in `ids`, in
the order given. Individuals may repeat.
"""
function select_individuals(hap::Haplotype, ids::AbstractVector{<:Integer})
    nid = hap.nhp ÷ 2
    all(i -> 1 <= i <= nid, ids) ||
        throw(ArgumentError("individuals must lie in 1:$nid"))
    isempty(ids) && throw(ArgumentError("select at least one individual"))
    cols = Vector{Int}(undef, 2length(ids))
    for (k, i) in enumerate(ids)
        cols[2k-1] = 2i - 1
        cols[2k] = 2i
    end
    return Haplotype(hap.gt[1:hap.nlc, cols])
end

"""
    select_loci(hap::Haplotype, loci) -> Haplotype

A new `Haplotype` with the rows `loci` (an index vector, a range or a
[`LocusSet`](@ref)), in the order given. Subset the [`VariantMap`](@ref)
with the same indices to keep the two aligned.
"""
function select_loci(hap::Haplotype, loci::AbstractVector{<:Integer})
    all(l -> 1 <= l <= hap.nlc, loci) ||
        throw(ArgumentError("loci must lie in 1:$(hap.nlc)"))
    isempty(loci) && throw(ArgumentError("select at least one locus"))
    return Haplotype(hap.gt[loci, :])
end
select_loci(hap::Haplotype, ls::LocusSet) = select_loci(hap, ls.loci)

"""
    dosage(g::Union{Haplotype,Genotype}, locus::Integer) -> Vector{Int8}
    dosage(g::Union{Haplotype,Genotype}, loci) -> Matrix{Int8}

Number of copies (0, 1 or 2) of the alternate allele carried by each
individual at `locus`, or at each of `loci` (an index vector, a range or a
[`LocusSet`](@ref)) as a `length(loci) × nid` matrix, one column per
individual.
"""
function dosage(g::Union{Haplotype,Genotype}, locus::Integer)
    1 <= locus <= _nlc(g) || throw(ArgumentError("locus must lie in 1:$(_nlc(g))"))
    return [Int8(_dosage(g, locus, ind)) for ind = 1:_nid(g)]
end

function dosage(g::Union{Haplotype,Genotype}, loci::AbstractVector{<:Integer})
    all(l -> 1 <= l <= _nlc(g), loci) || throw(ArgumentError("loci must lie in 1:$(_nlc(g))"))
    nid = _nid(g)
    Z = Matrix{Int8}(undef, length(loci), nid)
    Threads.@threads for ind = 1:nid
        @inbounds for (k, l) in enumerate(loci)
            Z[k, ind] = _dosage(g, l, ind)
        end
    end
    return Z
end
dosage(g::Union{Haplotype,Genotype}, ls::LocusSet) = dosage(g, ls.loci)
