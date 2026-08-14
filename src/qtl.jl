"""
    struct TraitQTL
Sparse representation of QTL loci and their effect sizes for a single trait.

# Fields
- `name::String`: Trait name.
- `loci::Vector{Int32}`: 1-based indices of QTL loci into a `VariantMap`.
- `additive::Vector{Float64}`: Additive effect sizes (α) per QTL locus.
- `dominance::Vector{Float64}`: Optional dominance effect sizes (d) per QTL locus (empty if purely additive).
"""
struct TraitQTL
    name::String
    loci::Vector{Int32}
    additive::Vector{Float64}
    dominance::Vector{Float64}

    function TraitQTL(
        name::AbstractString,
        loci::Vector{<:Integer},
        additive::Vector{<:Real};
        dominance::Vector{<:Real} = Float64[],
    )
        n = length(loci)
        length(additive) == n || throw(ArgumentError("loci and additive effect vectors must have the same length"))
        if !isempty(dominance)
            length(dominance) == n || throw(ArgumentError("dominance effect vector must match length of loci"))
        end
        all(loci .> 0) || throw(ArgumentError("Locus indices must be positive integers"))
        allunique(loci) || throw(ArgumentError("Locus indices must be unique"))
        new(String(name), Int32.(loci), Float64.(additive), Float64.(dominance))
    end
end

Base.length(qtl::TraitQTL) = length(qtl.loci)

function Base.show(io::IO, qtl::TraitQTL)
    print(io, "TraitQTL \"$(qtl.name)\" with $(length(qtl)) QTL loci")
end

"""
    struct MultiTraitQTL
Joint multi-trait QTL architecture containing a shared set of unique QTL loci
and a (K × T) effect matrix across T traits.

Enables fast BLAS matrix multiplications for calculating True Breeding Values (TBV)
across multiple traits simultaneously.
"""
struct MultiTraitQTL
    trait_names::Vector{String}
    loci::Vector{Int32}
    additive::Matrix{Float64}   # (n_qtl × n_traits)
    dominance::Matrix{Float64}  # (n_qtl × n_traits), optional (empty if none)

    function MultiTraitQTL(
        trait_names::Vector{<:AbstractString},
        loci::Vector{<:Integer},
        additive::Matrix{<:Real};
        dominance::Matrix{<:Real} = Matrix{Float64}(undef, 0, 0),
    )
        k = length(loci)
        t = length(trait_names)
        size(additive) == (k, t) ||
            throw(ArgumentError("additive effect matrix size $(size(additive)) must match (length(loci)=$k, length(traits)=$t)"))
        if !isempty(dominance)
            size(dominance) == (k, t) ||
                throw(ArgumentError("dominance effect matrix size must match (length(loci)=$k, length(traits)=$t)"))
        end
        all(loci .> 0) || throw(ArgumentError("Locus indices must be positive integers"))
        allunique(loci) || throw(ArgumentError("Locus indices must be unique"))
        new(String.(trait_names), Int32.(loci), Float64.(additive), Float64.(dominance))
    end
end

Base.length(mqtl::MultiTraitQTL) = length(mqtl.loci)

"""
    MultiTraitQTL(traits::Vector{TraitQTL})
Construct a unified `MultiTraitQTL` by merging multiple single-trait `TraitQTL` objects,
pooling all unique QTL loci into a shared set.
"""
function MultiTraitQTL(traits::Vector{TraitQTL})
    isempty(traits) && throw(ArgumentError("traits vector cannot be empty"))
    trait_names = [t.name for t in traits]
    
    # Collect all unique locus indices
    all_loci = sort(unique(vcat([t.loci for t in traits]...)))
    k = length(all_loci)
    t = length(traits)
    
    locus_to_row = Dict{Int32, Int}(all_loci[i] => i for i = 1:k)
    additive = zeros(Float64, k, t)
    
    has_dom = any(!isempty(t.dominance) for t in traits)
    dominance = has_dom ? zeros(Float64, k, t) : Matrix{Float64}(undef, 0, 0)
    
    for (trait_idx, tr) in enumerate(traits)
        for (i, locus) in enumerate(tr.loci)
            row = locus_to_row[locus]
            additive[row, trait_idx] = tr.additive[i]
            if has_dom && !isempty(tr.dominance)
                dominance[row, trait_idx] = tr.dominance[i]
            end
        end
    end
    
    return MultiTraitQTL(trait_names, all_loci, additive; dominance = dominance)
end

function Base.show(io::IO, mqtl::MultiTraitQTL)
    print(io, "MultiTraitQTL with $(length(mqtl.trait_names)) traits and $(length(mqtl)) unique QTL loci")
end

#
# True Breeding Value (TBV) Calculations
#

"""
    tbv(hap::Haplotype, qtl::TraitQTL) -> Vector{Float64}
Calculate True Breeding Values (TBV) for all individuals from a `Haplotype` (nlc × 2nid)
and a `TraitQTL`.
"""
function tbv(hap::Haplotype, qtl::TraitQTL)
    nid = hap.nhp ÷ 2
    tbvs = zeros(Float64, nid)
    loci = qtl.loci
    eff = qtl.additive
    dom = qtl.dominance
    has_dom = !isempty(dom)
    k = length(loci)

    @inbounds for ind = 1:nid
        h1_col = 2ind - 1
        h2_col = 2ind
        val = 0.0
        for j = 1:k
            loc = loci[j]
            a1 = Int(hap.gt[loc, h1_col])
            a2 = Int(hap.gt[loc, h2_col])
            dosage = a1 + a2
            val += dosage * eff[j]
            if has_dom && (dosage == 1)
                val += dom[j]
            end
        end
        tbvs[ind] = val
    end
    return tbvs
end

"""
    tbv(gt::Genotype, qtl::TraitQTL) -> Vector{Float64}
Calculate True Breeding Values (TBV) for all individuals from a `Genotype` (nid × 2nlc)
and a `TraitQTL`.
"""
function tbv(gt::Genotype, qtl::TraitQTL)
    nid = gt.nid
    tbvs = zeros(Float64, nid)
    loci = qtl.loci
    eff = qtl.additive
    dom = qtl.dominance
    has_dom = !isempty(dom)
    k = length(loci)

    @inbounds for ind = 1:nid
        val = 0.0
        for j = 1:k
            loc = loci[j]
            a1 = Int(gt.gt[ind, 2loc - 1])
            a2 = Int(gt.gt[ind, 2loc])
            dosage = a1 + a2
            val += dosage * eff[j]
            if has_dom && (dosage == 1)
                val += dom[j]
            end
        end
        tbvs[ind] = val
    end
    return tbvs
end

"""
    tbv(hap::Haplotype, mqtl::MultiTraitQTL) -> Matrix{Float64}
Calculate True Breeding Values (TBV) for all individuals across multiple traits (nid × n_traits)
from a `Haplotype` and a `MultiTraitQTL`.
"""
function tbv(hap::Haplotype, mqtl::MultiTraitQTL)
    nid = hap.nhp ÷ 2
    n_traits = length(mqtl.trait_names)
    k = length(mqtl.loci)
    
    # Extract dosage submatrix (nid × k)
    Z = Matrix{Float64}(undef, nid, k)
    loci = mqtl.loci
    
    @inbounds for j = 1:k
        loc = loci[j]
        for ind = 1:nid
            Z[ind, j] = Float64(Int(hap.gt[loc, 2ind - 1]) + Int(hap.gt[loc, 2ind]))
        end
    end
    
    # Linear algebra BLAS multiplication: (nid × k) * (k × n_traits) -> (nid × n_traits)
    tbvs = Z * mqtl.additive
    
    # Add dominance if present
    if !isempty(mqtl.dominance)
        @inbounds for ind = 1:nid
            for j = 1:k
                if Z[ind, j] == 1.0
                    for tr = 1:n_traits
                        tbvs[ind, tr] += mqtl.dominance[j, tr]
                    end
                end
            end
        end
    end
    
    return tbvs
end

"""
    tbv(gt::Genotype, mqtl::MultiTraitQTL) -> Matrix{Float64}
Calculate True Breeding Values (TBV) for all individuals across multiple traits (nid × n_traits)
from a `Genotype` and a `MultiTraitQTL`.
"""
function tbv(gt::Genotype, mqtl::MultiTraitQTL)
    nid = gt.nid
    n_traits = length(mqtl.trait_names)
    k = length(mqtl.loci)
    
    Z = Matrix{Float64}(undef, nid, k)
    loci = mqtl.loci
    
    @inbounds for j = 1:k
        loc = loci[j]
        for ind = 1:nid
            Z[ind, j] = Float64(Int(gt.gt[ind, 2loc - 1]) + Int(gt.gt[ind, 2loc]))
        end
    end
    
    tbvs = Z * mqtl.additive
    
    if !isempty(mqtl.dominance)
        @inbounds for ind = 1:nid
            for j = 1:k
                if Z[ind, j] == 1.0
                    for tr = 1:n_traits
                        tbvs[ind, tr] += mqtl.dominance[j, tr]
                    end
                end
            end
        end
    end
    
    return tbvs
end
