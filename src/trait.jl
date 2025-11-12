using Distributions

"""
    AbstractTrait
Abstract type for all traits.
"""
abstract type AbstractTrait end

"""
    abstract type Trait <: AbstractTrait end # non-threshold traits
Non-threshold traits have a field μ. 
"""
abstract type Trait <: AbstractTrait end # non-threshold traits

"""
    struct adTrait <: Trait
        name::AbstractString
    end
Trait with additive and dominance effects. To be implemented later.
"""
struct adTrait <: Trait
    name::AbstractString
end

"""
    struct adiTrait <: Trait
        name::AbstractString
    end
Trait with additive, dominance, and epistatic effects. To be implemented later.
"""
struct adiTrait <: Trait
    name::AbstractString
end

"""
    struct dtTrait <: AbstractTrait
        name::AbstractString
    end
Threshold trait with dominance effects. To be implemented later.
"""
struct dtTrait <: AbstractTrait
    name::AbstractString
end

"""
    struct ditTrait <: AbstractTrait
        name::AbstractString
    end
Threshold trait with dominance, and epistatic effects. To be implemented later.
"""
struct ditTrait <: AbstractTrait
    name::AbstractString
end

"""
    struct aTrait
        name::AbstractString # must be a valid Julia identifier
        sex::Int     # expressed in ♀(0), or ♂(1), or 2 for both sexes.
        age::Float64 # age that trait is measured
        h²::Float64  # narrow-sense heritability
        QTL::Symbol  # QTL column name in lmp
        μ::Float64   # population mean
        σₐ::Float64  # TBV std
        da::Distribution # QTL additive effect distribution, e.g., Normal()
        rev::Bool        # default true, higher the better
    end

Pure additive trait.
"""
struct aTrait <: Trait
    name::AbstractString # must be a valid Julia identifier
    sex::Int     # expressed in ♀(0), or ♂(1), or 2 for both sexes
    age::Float64 # age that trait is measured
    h²::Float64
    QTL::Symbol  # QTL column name in lmp
    μ::Float64   # population mean
    σₐ::Float64  # TBV std
    da::Distribution # QTL additive effect distribution, e.g., Normal()
    rev::Bool        # default true, higher the better
end

"""
    struct tTrait <: AbstractTrait
        name::AbstractString # must be a valid Julia identifier
        threshold::Vector{Float64} # phenotype will be 0, 1, ...
        sex::Int     # expressed in ♀(0), or ♂(1), 2 for both sexes.
        age::Float64 # age that trait is measured
        h²::Float64  # narrow-sense heritability
        QTL::Symbol
        σₐ::Float64  # TBV std
        da::Distribution  # e.g., Normal()
        rev::Bool    # default true, higher the better
    end
Struct for threshold traits. The number of phenotype categories equals the
length of threshold + 1. The underlying genotype is continuous. The phenotype
will be one of 0, 1, ..., `length(threshold)`. Population mean is not used for
threshold traits.
"""
struct tTrait <: AbstractTrait
    name::AbstractString # must be a valid Julia identifier
    threshold::Vector{Float64} # phenotype will be 0, 1, ...
    sex::Int     # expressed in ♀(0), or ♂(1), 2 for both sexes
    age::Float64 # age that trait is measured
    h²::Float64  # narrow-sense heritability
    QTL::Symbol
    σₐ::Float64  # TBV std
    da::Distribution  # e.g., Normal()
    rev::Bool    # default true, higher the better
end

"""
    Trait(
        name;
        sex = 2,
        age = 0.,
        h² = 0.25,
        QTL = :qtl,
        μ = 0.0,
        σₐ = 1.0,
        da = Normal(),
        rev = true,
    )
Helper function to create a Trait object. It checks the validity of the input
parameters and returns a Trait object. Refer to the `Trait` struct for the
description of the parameters.
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
    rev = true,
)
    occursin(r"^[[:alpha:]_\p{L}][\p{L}\p{N}_]*$", name) ||
        error("Name must be valid for a data frame column name")
    sex ∈ 0:2 || error("Sex can only be in 0:2")
    age ≥ 0 || error("Age must be non negative")
    0.0 < h² ≤ 1.0 || error("h² must be in (0, 1]")
    σₐ > 0 || error("σₐ must be positive")
    aTrait(name, sex, age, h², QTL, μ, σₐ, da, rev)
end

"""
    Trait(
        name::AbstractString,
        weight::Vector{Float64};
        sex = 2,
        age = 1.0,
        h² = 0.25,
        QTL = :qtl,
        σₐ = 1.0,
        da = Normal(),
        rev = true,
    )
Helper function to create a threshold `tTrait`. It checks the validity of the
input parameters and returns a `tTrait`. Refer to the `tTrait` struct for the
description of the parameters. The weight, of minimum length 2, defines the
relative ratio of the phenotype categories. The length of the weight is one more
than the length of the thresholds in struct `tTrait`. This function then decides
the thresholds based on the weight and standard normal distribution of the
genotype. 

The underlying TBV distribution is Normal(0, 1.0). The categorical phenotypes
are calculated from the phenotypic thresholds calculated from `Weights`.
"""
function Trait(
    name::AbstractString,
    weight::Vector{Float64};
    sex = 2,
    age = 1.0,
    h² = 0.25,
    QTL = :qtl,
    da = Normal(),
    rev = true,
)
    occursin(r"^[[:alpha:]_\p{L}][\p{L}\p{N}_]*$", name) ||
        error("Name must be valid for a data frame column name")
    all(weight .> 0) || error("Weights should be all > 0")
    sex ∈ 0:2 || error("Sex can only be in 0:2")
    age ≥ 0 || error("Age must be non negative")
    0.0 < h² ≤ 1.0 || error("h² must be in (0, 1]")
    weight ./= sum(weight)
    σₚ = sqrt(1.0 / h²)
    threshold = map(x -> quantile(Normal(0, σₚ), x), cumsum(weight[1:(end-1)]))
    tTrait(name, threshold, sex, age, h², QTL, 1.0, da, rev)
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
    isa(trt, Trait) && println(io, "  Population mean: $(trt.μ)")
    println(io, "     QTL set name: $(trt.QTL)")
    println(io, "         std(TBV): $(trt.σₐ)")
    println(io, " QTL add. distri.: $(trt.da)")
    print(io, "Higher the better: $(trt.rev)")
end
