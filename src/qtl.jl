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

# Constructors
    TraitQTL(name::AbstractString, loci::Vector{<:Integer}, additive::Vector{<:Real}; dominance::Vector{<:Real} = Float64[]) -> TraitQTL

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
and effect matrices across multiple traits.

Enables high-performance BLAS matrix multiplications for calculating True Breeding Values
(TBV) across all traits and individuals simultaneously.

# Fields
- `trait_names::Vector{String}`: Names of the traits.
- `loci::Vector{Int32}`: Sorted unique 1-based indices of all union QTL loci into a [`VariantMap`](@ref).
- `additive::Matrix{Float64}`: Additive effect matrix of dimension `(length(loci), length(trait_names))`.
- `dominance::Matrix{Float64}`: Dominance effect matrix of dimension `(length(loci), length(trait_names))` (or empty `0 × 0` matrix).

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

Calculate True Breeding Values (TBV) for all individuals from a locus-major `Haplotype`
matrix and a single-trait [`TraitQTL`](@ref).

For each individual `i`, dosage `Z_j = h[j, 2i-1] + h[j, 2i] ∈ {0, 1, 2}` is evaluated
at each QTL locus `j`. The individual's genetic value is:
```math
\\text{TBV}_i = \\sum_{j=1}^K Z_{ij} \\alpha_j + \\sum_{j: Z_{ij} = 1} d_j
```
where `α` denotes additive effect sizes and `d` denotes dominance deviations (for heterozygous loci).

# Arguments
- `hap::Haplotype`: Phased haplotype matrix of dimension `(nlc, 2 * nid)`.
- `qtl::TraitQTL`: Single-trait QTL architecture.

# Returns
- `Vector{Float64}`: Vector of length `nid` containing TBVs for each individual.

# Examples
```julia
using BnGStructs

hps = Haplotype(10, 4) # 10 loci, 2 individuals
hps[1, 1] = true; hps[1, 2] = true # individual 1 homozygous at locus 1 (dosage = 2)
qtl = TraitQTL("Trait1", [1], [1.5])
tbv(hps, qtl) # [3.0, 0.0]
```
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

Calculate True Breeding Values (TBV) for all individuals from an individual-major `Genotype`
matrix and a single-trait [`TraitQTL`](@ref).

For each individual `i` and QTL locus `l_j`, allele dosage `Z_ij = g[i, 2l_j-1] + g[i, 2l_j] ∈ {0, 1, 2}`.
The individual's genetic value is:
```math
\\text{TBV}_i = \\sum_{j=1}^K Z_{ij} \\alpha_j + \\sum_{j: Z_{ij} = 1} d_j
```

# Arguments
- `gt::Genotype`: Phased genotype matrix of dimension `(nid, 2 * nlc)`.
- `qtl::TraitQTL`: Single-trait QTL architecture.

# Returns
- `Vector{Float64}`: Vector of length `nid` containing TBVs for each individual.

# Examples
```julia
using BnGStructs

gt = Genotype(2, 20) # 2 individuals, 10 loci
gt[1, 1] = true; gt[1, 2] = true # individual 1 homozygous at locus 1 (dosage = 2)
qtl = TraitQTL("Trait1", [1], [1.5])
tbv(gt, qtl) # [3.0, 0.0]
```
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

Calculate True Breeding Values (TBV) for all individuals across multiple traits from
a locus-major `Haplotype` matrix and a [`MultiTraitQTL`](@ref).

Extracts the dosage submatrix `Z` of size `(nid, K)` across the `K` union QTL loci and
computes the additive component via optimized BLAS matrix multiplication:
```math
\\mathbf{TBV} = Z \\mathbf{A} + \\mathbf{D}
```
where `A` is the `(K × T)` additive effect matrix and `D` accounts for dominance
deviations on heterozygous genotypes (`Z_ij = 1`).

# Arguments
- `hap::Haplotype`: Phased haplotype matrix of dimension `(nlc, 2 * nid)`.
- `mqtl::MultiTraitQTL`: Joint multi-trait QTL architecture across `T` traits.

# Returns
- `Matrix{Float64}`: Matrix of dimension `(nid, T)` containing TBVs for all individuals across all traits.
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

Calculate True Breeding Values (TBV) for all individuals across multiple traits from
an individual-major `Genotype` matrix and a [`MultiTraitQTL`](@ref).

Extracts the dosage submatrix `Z` of size `(nid, K)` across the `K` union QTL loci and
computes the additive component via optimized BLAS matrix multiplication:
```math
\\mathbf{TBV} = Z \\mathbf{A} + \\mathbf{D}
```
where `A` is the `(K × T)` additive effect matrix and `D` accounts for dominance
deviations on heterozygous genotypes (`Z_ij = 1`).

# Arguments
- `gt::Genotype`: Phased genotype matrix of dimension `(nid, 2 * nlc)`.
- `mqtl::MultiTraitQTL`: Joint multi-trait QTL architecture across `T` traits.

# Returns
- `Matrix{Float64}`: Matrix of dimension `(nid, T)` containing TBVs for all individuals across all traits.
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
