"""
    struct TraitQTL

Sparse representation of quantitative trait loci (QTL) and genetic effect sizes for a single trait.

Causal loci are identified by 1-based indices into a corresponding [`VariantMap`](@ref).
Genetic effects at each locus can include additive allele substitution effects (α) and optional
dominance deviations (d).

# Fields
- `name::String`: Trait name.
- `loci::Vector{Int32}`: Unique 1-based indices of QTL loci into a [`VariantMap`](@ref).
- `additive::Vector{Float64}`: Additive effect sizes (α) per QTL locus.
- `dominance::Vector{Float64}`: Optional dominance effect sizes (d) per QTL locus (empty if purely additive).
- `intercept::Float64`: Constant added to every genotypic value (default `0.0`). Use it to
  place the base-population mean without shifting the sampled effects.

Genotype values at a locus are `0`, `a + d` and `2a` for 0, 1 and 2 copies
of the alternate allele, so `a` is half the difference between the
homozygotes and `d` the deviation of the heterozygote from their midpoint.

# Constructors
    TraitQTL(name, loci, additive; dominance = Float64[], intercept = 0.0) -> TraitQTL

Construct a single-trait QTL architecture. Requires that `loci` and `additive` vectors have
the same length, all locus indices are positive and unique, and that `dominance` (if provided)
matches the length of `loci`.

# Methods
- `length(qtl::TraitQTL)`: Number of QTL loci for the trait.

# Examples
```julia
using BnGStructs

# Pure additive trait with 2 QTLs
qtl1 = TraitQTL("Yield", [10, 45], [0.5, -0.3])
length(qtl1) # 2

# Trait with additive and dominance effects
qtl2 = TraitQTL("Growth", [10, 20], [1.2, 0.8]; dominance=[0.4, -0.1])
```
"""
struct TraitQTL
    name::String
    loci::Vector{Int32}
    additive::Vector{Float64}
    dominance::Vector{Float64}
    intercept::Float64

    function TraitQTL(
        name::AbstractString,
        loci::AbstractVector{<:Integer},
        additive::AbstractVector{<:Real};
        dominance::AbstractVector{<:Real} = Float64[],
        intercept::Real = 0.0,
    )
        n = length(loci)
        length(additive) == n || throw(ArgumentError("loci and additive effect vectors must have the same length"))
        if !isempty(dominance)
            length(dominance) == n || throw(ArgumentError("dominance effect vector must match length of loci"))
        end
        all(loci .> 0) || throw(ArgumentError("Locus indices must be positive integers"))
        allunique(loci) || throw(ArgumentError("Locus indices must be unique"))
        isfinite(intercept) || throw(ArgumentError("intercept must be finite"))
        new(String(name), Int32.(loci), Float64.(additive), Float64.(dominance), Float64(intercept))
    end
end

Base.length(qtl::TraitQTL) = length(qtl.loci)

function Base.show(io::IO, qtl::TraitQTL)
    print(io, "TraitQTL \"$(qtl.name)\" with $(length(qtl)) QTL loci")
end

"""
    struct MultiTraitQTL

Joint multi-trait QTL architecture containing a shared set of unique QTL loci
and effect matrices across multiple traits.

Enables high-performance BLAS matrix multiplications for calculating True Breeding Values
(TBV) across all traits and individuals simultaneously.

# Fields
- `trait_names::Vector{String}`: Names of the traits.
- `loci::Vector{Int32}`: Sorted unique 1-based indices of all union QTL loci into a [`VariantMap`](@ref).
- `additive::Matrix{Float64}`: Additive effect matrix of dimension `(length(loci), length(trait_names))`.
- `dominance::Matrix{Float64}`: Dominance effect matrix of dimension `(length(loci), length(trait_names))` (or empty `0 × 0` matrix).
- `intercept::Vector{Float64}`: One constant per trait, added to every genotypic value (default zeros).

# Constructors
- `MultiTraitQTL(trait_names::Vector{<:AbstractString}, loci::Vector{<:Integer}, additive::Matrix{<:Real}; dominance = Matrix{Float64}(undef, 0, 0))`
- `MultiTraitQTL(traits::Vector{TraitQTL})`: Merges multiple single-trait `TraitQTL` objects into a unified architecture.

# Examples
```julia
using BnGStructs

qtl1 = TraitQTL("Trait1", [1, 2], [0.5, -1.0])
qtl2 = TraitQTL("Trait2", [2, 3], [2.0, 3.0])
mqtl = MultiTraitQTL([qtl1, qtl2])

length(mqtl.trait_names) # 2
length(mqtl)             # 3 unique QTL loci
size(mqtl.additive)      # (3, 2)
```
"""
struct MultiTraitQTL
    trait_names::Vector{String}
    loci::Vector{Int32}
    additive::Matrix{Float64}   # (n_qtl × n_traits)
    dominance::Matrix{Float64}  # (n_qtl × n_traits), optional (empty if none)
    intercept::Vector{Float64}  # one per trait

    function MultiTraitQTL(
        trait_names::AbstractVector{<:AbstractString},
        loci::AbstractVector{<:Integer},
        additive::AbstractMatrix{<:Real};
        dominance::AbstractMatrix{<:Real} = Matrix{Float64}(undef, 0, 0),
        intercept::AbstractVector{<:Real} = zeros(length(trait_names)),
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
        length(intercept) == t && all(isfinite, intercept) ||
            throw(ArgumentError("intercept needs one finite value per trait"))
        new(String.(trait_names), Int32.(loci), Float64.(additive), Float64.(dominance), Float64.(intercept))
    end
end

Base.length(mqtl::MultiTraitQTL) = length(mqtl.loci)

"""
    MultiTraitQTL(traits::Vector{TraitQTL}) -> MultiTraitQTL

Construct a unified `MultiTraitQTL` by merging a vector of single-trait [`TraitQTL`](@ref) objects.

The union of all unique QTL locus indices across traits is identified and sorted. Additive and
dominance effect matrices are assembled such that row `k` corresponds to the `k`-th unique locus and
column `t` corresponds to the `t`-th trait. Any locus not present in a given trait is assigned an effect of `0.0`.

# Arguments
- `traits::Vector{TraitQTL}`: Non-empty list of single-trait QTL definitions.

# Returns
- `MultiTraitQTL`: Combined multi-trait QTL architecture.
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
    
    return MultiTraitQTL(
        trait_names, all_loci, additive;
        dominance = dominance, intercept = [t.intercept for t in traits],
    )
end

function Base.show(io::IO, mqtl::MultiTraitQTL)
    print(io, "MultiTraitQTL with $(length(mqtl.trait_names)) traits and $(length(mqtl)) unique QTL loci")
end

#
# Genotypic values, additive genetic values and breeding values
#

# Dosage of the alternate allele of individual `ind` at locus `loc`.
@inline _dosage(hap::Haplotype, loc, ind) =
    @inbounds Int(hap.gt[loc, 2ind-1]) + Int(hap.gt[loc, 2ind])
@inline _dosage(gt::Genotype, loc, ind) =
    @inbounds Int(gt.gt[ind, 2loc-1]) + Int(gt.gt[ind, 2loc])
_nid(hap::Haplotype) = hap.nhp ÷ 2
_nid(gt::Genotype) = gt.nid
_nlc(hap::Haplotype) = hap.nlc
_nlc(gt::Genotype) = gt.nas ÷ 2

function _check_loci(g::Union{Haplotype,Genotype}, loci)
    isempty(loci) || maximum(loci) <= _nlc(g) ||
        throw(ArgumentError("QTL locus $(maximum(loci)) is beyond the $(_nlc(g)) loci of the genotypes"))
    return nothing
end

"""
    genotypic_value(g::Union{Haplotype,Genotype}, qtl::TraitQTL) -> Vector{Float64}
    genotypic_value(g::Union{Haplotype,Genotype}, qtl::MultiTraitQTL) -> Matrix{Float64}

Genotypic value of every individual: with dosage ``z_{ij} \\in \\{0, 1, 2\\}``
of the alternate allele of individual ``i`` at QTL ``j``,

```math
G_i = \\mu + \\sum_j z_{ij} a_j + \\sum_{j:\\, z_{ij} = 1} d_j ,
```

where ``\\mu`` is `qtl.intercept`, ``a`` is `qtl.additive` and ``d`` is
`qtl.dominance`. Loci are independent in their effects (no epistasis). The
multi-trait method returns an `nid × ntraits` matrix.

With no dominance this is the additive genetic value, which [`tbv`](@ref)
returns under its usual name. With dominance it is **not** a breeding
value; see [`breeding_value`](@ref).

# Examples
```julia
using BnGStructs

hps = Haplotype(10, 4)                  # 10 loci, 2 individuals
hps[1, 1] = true; hps[1, 2] = true      # individual 1 is homozygous at locus 1
hps[1, 3] = true                        # individual 2 is heterozygous
qtl = TraitQTL("T", [1], [1.5]; dominance = [0.5])
genotypic_value(hps, qtl)               # [3.0, 2.0]
```
"""
function genotypic_value(g::Union{Haplotype,Genotype}, qtl::TraitQTL)
    loci, eff, dom = qtl.loci, qtl.additive, qtl.dominance
    _check_loci(g, loci)
    has_dom = !isempty(dom)
    vals = fill(qtl.intercept, _nid(g))
    @inbounds for ind in eachindex(vals)
        val = 0.0
        for j in eachindex(loci)
            z = _dosage(g, loci[j], ind)
            val += z * eff[j]
            has_dom && z == 1 && (val += dom[j])
        end
        vals[ind] += val
    end
    return vals
end

function genotypic_value(g::Union{Haplotype,Genotype}, mqtl::MultiTraitQTL)
    loci = mqtl.loci
    _check_loci(g, loci)
    nid, k = _nid(g), length(loci)
    Z = Matrix{Float64}(undef, nid, k)
    @inbounds for j = 1:k, ind = 1:nid
        Z[ind, j] = _dosage(g, loci[j], ind)
    end
    vals = Z * mqtl.additive
    if !isempty(mqtl.dominance)
        vals .+= (Z .== 1.0) * mqtl.dominance
    end
    vals .+= mqtl.intercept'
    return vals
end

_has_dominance(qtl::TraitQTL) = any(!iszero, qtl.dominance)
_has_dominance(qtl::MultiTraitQTL) = any(!iszero, qtl.dominance)

"""
    tbv(g::Union{Haplotype,Genotype}, qtl::TraitQTL) -> Vector{Float64}
    tbv(g::Union{Haplotype,Genotype}, qtl::MultiTraitQTL) -> Matrix{Float64}

True breeding value under an **additive** architecture:
``\\mu + \\sum_j z_{ij} a_j``, the same as [`genotypic_value`](@ref). The
values are not centred; they include `qtl.intercept` and the mean of the
population at hand.

Throws an `ArgumentError` if `qtl` has non-zero dominance effects, because
the sum above is then neither the genotypic value nor the breeding value.
Use [`genotypic_value`](@ref) or [`breeding_value`](@ref) instead.

# Examples
```julia
using BnGStructs

hps = Haplotype(10, 4) # 10 loci, 2 individuals
hps[1, 1] = true; hps[1, 2] = true # individual 1 homozygous at locus 1 (dosage = 2)
qtl = TraitQTL("Trait1", [1], [1.5])
tbv(hps, qtl) # [3.0, 0.0]
```
"""
function tbv(g::Union{Haplotype,Genotype}, qtl::Union{TraitQTL,MultiTraitQTL})
    _has_dominance(qtl) && throw(ArgumentError(
        "tbv is defined for additive architectures; with dominance use " *
        "genotypic_value (G) or breeding_value (A)"))
    return genotypic_value(g, qtl)
end

"""
    breeding_value(g, qtl::TraitQTL; p = nothing) -> Vector{Float64}
    breeding_value(g, qtl::MultiTraitQTL; p = nothing) -> Matrix{Float64}

Breeding value (additive genetic value) of every individual, as a
deviation from the mean of the reference population:

```math
A_i = \\sum_j (z_{ij} - 2p_j)\\,\\alpha_j, \\qquad
\\alpha_j = a_j + d_j (1 - 2p_j),
```

where ``p_j`` is the frequency of the alternate allele at QTL ``j`` in the
reference population and ``\\alpha_j`` is the average effect of an allele
substitution there. `p` gives those frequencies, one per QTL in the order
of `qtl.loci`; by default they are the frequencies among the individuals
of `g`, which makes that sample the reference population.

Assumptions: a randomly mating reference population in Hardy–Weinberg
proportions, and effects that add across loci. Under linkage
disequilibrium the values are still the sum of single-locus breeding
values; their variance then includes covariances between loci.

With no dominance, `breeding_value(g, qtl)` equals
`tbv(g, qtl) .- mean(tbv(g, qtl))` when `p` is taken from `g`.
"""
function breeding_value(
    g::Union{Haplotype,Genotype},
    qtl::TraitQTL;
    p::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    loci = qtl.loci
    _check_loci(g, loci)
    nid = _nid(g)
    freq = _qtl_frequencies(g, loci, p)
    vals = zeros(nid)
    has_dom = !isempty(qtl.dominance)
    @inbounds for j in eachindex(loci)
        pj = freq[j]
        α = qtl.additive[j] + (has_dom ? qtl.dominance[j] * (1 - 2pj) : 0.0)
        for ind = 1:nid
            vals[ind] += (_dosage(g, loci[j], ind) - 2pj) * α
        end
    end
    return vals
end

function breeding_value(
    g::Union{Haplotype,Genotype},
    mqtl::MultiTraitQTL;
    p::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    loci = mqtl.loci
    _check_loci(g, loci)
    nid, k = _nid(g), length(loci)
    freq = _qtl_frequencies(g, loci, p)
    Z = Matrix{Float64}(undef, nid, k)
    @inbounds for j = 1:k, ind = 1:nid
        Z[ind, j] = _dosage(g, loci[j], ind) - 2freq[j]
    end
    α = isempty(mqtl.dominance) ? mqtl.additive :
        mqtl.additive .+ mqtl.dominance .* (1 .- 2 .* freq)
    return Z * α
end

function _qtl_frequencies(g::Union{Haplotype,Genotype}, loci, p)
    if p !== nothing
        length(p) == length(loci) ||
            throw(DimensionMismatch("p needs one frequency per QTL ($(length(loci)))"))
        all(x -> 0 <= x <= 1, p) || throw(ArgumentError("frequencies must lie in [0, 1]"))
        return Float64.(p)
    end
    nid = _nid(g)
    return [sum(_dosage(g, l, ind) for ind = 1:nid) / 2nid for l in loci]
end
