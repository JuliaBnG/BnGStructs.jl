using Distributions

"""
    abstract type AbstractTrait

Abstract supertype for all genetic traits (quantitative and threshold/categorical).

# Subtypes
- [`Trait`](@ref): Abstract supertype for continuous traits with a population mean `μ` (e.g. [`aTrait`](@ref)).
- [`tTrait`](@ref): Concrete type for threshold/ordinal traits modeled via liability thresholds.
- Planned trait architectures: [`adTrait`](@ref), [`adiTrait`](@ref), [`dtTrait`](@ref), [`ditTrait`](@ref).
"""
abstract type AbstractTrait end

"""
    abstract type Trait <: AbstractTrait

Abstract supertype for continuous (non-threshold) quantitative traits.

Subtypes of `Trait` possess a population phenotypic mean field `μ`.
"""
abstract type Trait <: AbstractTrait end # non-threshold traits

"""
    struct adTrait <: Trait

Planned trait type with additive and dominance genetic effects.
"""
struct adTrait <: Trait
    name::String
end

"""
    struct adiTrait <: Trait

Planned trait type with additive, dominance, and epistatic interaction effects.
"""
struct adiTrait <: Trait
    name::String
end

"""
    struct dtTrait <: AbstractTrait

Planned threshold trait type with dominance genetic effects.
"""
struct dtTrait <: AbstractTrait
    name::String
end

"""
    struct ditTrait <: AbstractTrait

Planned threshold trait type with dominance and epistatic interaction effects.
"""
struct ditTrait <: AbstractTrait
    name::String
end

"""
    struct aTrait{D<:Distribution} <: Trait

A continuous quantitative trait with purely additive genetic architecture.

# Fields
- `name::String`: Trait name (must be a valid Julia identifier).
- `sex::Int`: Sex-limited expression: `0` for females (♀), `1` for males (♂), or `2` for both sexes.
- `age::Float64`: Measurement age or developmental stage (non-negative).
- `h²::Float64`: Narrow-sense heritability (`0 < h² ≤ 1`).
- `QTL::Symbol`: Column name identifying corresponding QTL loci in a locus map.
- `μ::Float64`: Population phenotypic mean.
- `σₐ::Float64`: Additive genetic standard deviation (standard deviation of true breeding values, `σₐ > 0`).
- `da::D`: Probability distribution used to sample QTL additive effect sizes (e.g., `Distributions.Normal()`).

See also: [`Trait`](@ref), [`tTrait`](@ref).
"""
struct aTrait{D<:Distribution} <: Trait
    name::String # must be a valid Julia identifier
    sex::Int     # expressed in ♀(0), or ♂(1), or 2 for both sexes
    age::Float64 # age that trait is measured
    h²::Float64
    QTL::Symbol  # QTL column name in lmp
    μ::Float64   # population mean
    σₐ::Float64  # TBV std
    da::D        # QTL additive effect distribution, e.g., Normal()
end

"""
    struct tTrait{D<:Distribution} <: AbstractTrait

A categorical or ordinal threshold trait modeled under the polygenic liability-threshold framework.

In threshold traits, the underlying liability follows a continuous additive distribution with standard deviation
`σₐ = 1.0`. Phenotypic categories (`0, 1, ..., K`) are demarcated by `K` real-valued thresholds,
where `K = length(threshold)`.

# Fields
- `name::String`: Trait name (must be a valid Julia identifier).
- `threshold::Vector{Float64}`: Sorted boundary cutoffs on the phenotypic liability scale.
- `sex::Int`: Sex-limited expression: `0` for females (♀), `1` for males (♂), or `2` for both sexes.
- `age::Float64`: Measurement age or developmental stage (non-negative).
- `h²::Float64`: Narrow-sense heritability on the liability scale (`0 < h² ≤ 1`).
- `QTL::Symbol`: Column name identifying corresponding QTL loci in a locus map.
- `σₐ::Float64`: Additive genetic standard deviation on the liability scale (fixed at `1.0`).
- `da::D`: Probability distribution used to sample QTL additive effect sizes on the liability scale.

See also: [`Trait`](@ref), [`aTrait`](@ref).
"""
struct tTrait{D<:Distribution} <: AbstractTrait
    name::String # must be a valid Julia identifier
    threshold::Vector{Float64} # phenotype will be 0, 1, ...
    sex::Int     # expressed in ♀(0), or ♂(1), 2 for both sexes
    age::Float64 # age that trait is measured
    h²::Float64  # narrow-sense heritability
    QTL::Symbol
    σₐ::Float64  # TBV std
    da::D        # e.g., Normal()
end

"""
    Trait(
        name::AbstractString;
        sex = 2,
        age = 0.0,
        h² = 0.25,
        QTL = :qtl,
        μ = 0.0,
        σₐ = 1.0,
        da = Normal(),
    ) -> aTrait

Construct a continuous additive quantitative trait ([`aTrait`](@ref)).

# Arguments
- `name::AbstractString`: Name of the trait. Must be a valid Julia identifier suitable for DataFrame column names.

# Keywords
- `sex`: Sex-limited expression: `0` for females (♀), `1` for males (♂), `2` for both (default: `2`).
- `age`: Measurement age or stage (default: `0.0`, must be finite and `>= 0`).
- `h²`: Narrow-sense heritability (default: `0.25`, must be in `(0, 1]`).
- `QTL`: Symbol naming the QTL marker column in locus maps (default: `:qtl`).
- `μ`: Population phenotypic mean (default: `0.0`, must be finite).
- `σₐ`: Additive genetic standard deviation (default: `1.0`, must be finite and `> 0`).
- `da`: Distribution of QTL additive effects (default: `Normal()`).

# Returns
- `aTrait`: Validated continuous trait definition.

# Examples
```julia
using BnGStructs, Distributions

milk = Trait("MilkYield"; h²=0.35, μ=30.0, σₐ=2.5, da=Normal(0, 1.2))
```
"""
function Trait(
    name::AbstractString;
    sex = 2,
    age = 0.0,
    h² = 0.25,
    QTL = :qtl,
    μ = 0.0,
    σₐ = 1.0,
    da = Normal(),
)
    occursin(r"^[[:alpha:]_\p{L}][\p{L}\p{N}_]*$", name) ||
        error("Name must be valid for a data frame column name")
    sex ∈ 0:2 || error("Sex can only be in 0:2")
    isfinite(age) && age ≥ 0 || error("Age must be finite and non-negative")
    isfinite(h²) && 0.0 < h² ≤ 1.0 || error("h² must be finite and in (0, 1]")
    isfinite(μ) || error("μ must be finite")
    isfinite(σₐ) && σₐ > 0 || error("σₐ must be finite and positive")
    aTrait(String(name), sex, Float64(age), Float64(h²), Symbol(QTL), Float64(μ), Float64(σₐ), da)
end

"""
    Trait(
        name::AbstractString,
        weight::AbstractVector{<:Real};
        sex = 2,
        age = 1.0,
        h² = 0.25,
        QTL = :qtl,
        da = Normal(),
    ) -> tTrait

Construct a threshold trait ([`tTrait`](@ref)) under the liability-threshold model.

The relative frequencies or proportions of the categorical phenotypic outcomes are specified
by `weight`. For `K + 1` phenotype categories (`0, 1, ..., K`), `weight` must contain at least 2
positive values. Phenotypic thresholds are computed via quantiles of `Normal(0, σₚ)`
where phenotypic standard deviation `σₚ = sqrt(1 / h²)`.

# Arguments
- `name::AbstractString`: Name of the trait (must be a valid Julia identifier).
- `weight::AbstractVector{<:Real}`: Relative proportions or weights for each category (length `>= 2`).

# Keywords
- `sex`: Sex-limited expression: `0` for females (♀), `1` for males (♂), `2` for both (default: `2`).
- `age`: Measurement age or stage (default: `1.0`, must be finite and `>= 0`).
- `h²`: Narrow-sense heritability on the liability scale (default: `0.25`, must be in `(0, 1]`).
- `QTL`: Symbol naming the QTL marker column in locus maps (default: `:qtl`).
- `da`: Distribution of QTL additive effects on the liability scale (default: `Normal()`).

# Returns
- `tTrait`: Validated threshold trait definition with computed thresholds.

# Examples
```julia
using BnGStructs

# Binary disease trait: 90% unaffected (0), 10% affected (1)
disease = Trait("Mastitis", [0.9, 0.1]; h²=0.10)
length(disease.threshold) # 1

# Three categories: low (50%), medium (30%), high (20%)
calving_ease = Trait("CalvingEase", [0.5, 0.3, 0.2]; h²=0.20)
length(calving_ease.threshold) # 2
```
"""
function Trait(
    name::AbstractString,
    weight::AbstractVector{<:Real};
    sex = 2,
    age = 1.0,
    h² = 0.25,
    QTL = :qtl,
    da = Normal(),
)
    occursin(r"^[[:alpha:]_\p{L}][\p{L}\p{N}_]*$", name) ||
        error("Name must be valid for a data frame column name")
    length(weight) ≥ 2 || error("At least two category weights are required")
    all(isfinite, weight) && all(>(0), weight) ||
        error("Weights must be finite and positive")
    sex ∈ 0:2 || error("Sex can only be in 0:2")
    isfinite(age) && age ≥ 0 || error("Age must be finite and non-negative")
    isfinite(h²) && 0.0 < h² ≤ 1.0 || error("h² must be finite and in (0, 1]")
    normalized_weight = Float64.(weight)
    total_weight = sum(normalized_weight)
    isfinite(total_weight) || error("Weights must have a finite sum")
    normalized_weight ./= total_weight
    σₚ = sqrt(1.0 / h²)
    threshold = map(
        x -> quantile(Normal(0, σₚ), x),
        cumsum(normalized_weight[1:(end-1)]),
    )
    tTrait(String(name), threshold, sex, Float64(age), Float64(h²), Symbol(QTL), 1.0, da)
end

function Base.show(io::IO, trt::AbstractTrait)
    println(io, "             Name: $(trt.name)")
    isa(trt, tTrait) && println(io, "       Thresholds: $(trt.threshold)")
    if trt.sex == 0
        println(io, "     Expresses in: females")
    elseif trt.sex == 1
        println(io, "       Express in: males")
    elseif trt.sex == 2
        println(io, "       Express in: both sexes")
    else
        println(io, "       Not expressed/measured")
    end
    println(io, "      Express age: $(trt.age)")
    println(io, "     Heritability: $(trt.h²)")
    isa(trt, Trait) && println(io, "  Init. pop. mean: $(trt.μ)")
    println(io, "     QTL set name: $(trt.QTL)")
    println(io, "         std(TBV): $(trt.σₐ)")
    println(io, " QTL add. distri.: $(trt.da)")
end
