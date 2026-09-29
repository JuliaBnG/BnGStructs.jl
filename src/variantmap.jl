"""
    struct VariantMap

A compact struct-of-arrays representation of genomic variant coordinates and alleles.

`VariantMap` stores per-variant metadata across parallel contiguous vectors for cache-efficient
access and minimal memory overhead.

# Fields
- `chr::Vector{Int8}`: Chromosome numbers (1-based integer identifiers).
- `pos::Vector{UInt32}`: Physical base-pair coordinates on the chromosome.
- `ref::Vector{Char}`: Reference alleles (e.g. `'A'`, `'C'`, `'G'`, `'T'`).
- `alt::Vector{Char}`: Alternate alleles.
- `frq::Vector{Float32}`: Alternate allele frequencies (optional; empty vector `Float32[]` if omitted).

# Constructors
- `VariantMap(chr::Vector{Int8}, pos::Vector{UInt32}, ref::Vector{Char}, alt::Vector{Char}, frq::Vector{Float32} = Float32[])`
- `VariantMap(chr::Vector{<:Integer}, pos::Vector{<:Integer}, ref::Vector{Char}, alt::Vector{Char}; frq = Float32[])`

# Indexing & Methods
- `length(vm)`: Returns the number of variants.
- `size(vm)`: Returns `(length(vm),)`.
- `vm[i]::NamedTuple`: Returns variant details `(chr, pos, ref, alt, frq)` at index `i`.
  If frequencies are omitted, `frq` is `NaN32`.
- `vm[idxs]::VariantMap`: Returns a new `VariantMap` containing the subset of variants at indices `idxs`.

# Examples
```julia
using BnGStructs

vm = VariantMap([1, 1, 2], [1000, 2000, 1500], ['A', 'C', 'G'], ['G', 'T', 'A']; frq=[0.1, 0.4, 0.25])
length(vm)       # 3
vm[1]            # (chr = 1, pos = 0x000003e8, ref = 'A', alt = 'G', frq = 0.1f0)
sub = vm[[1, 3]] # VariantMap with 2 variants
```
"""
struct VariantMap
    chr::Vector{Int8}
    pos::Vector{UInt32}
    ref::Vector{Char}
    alt::Vector{Char}
    frq::Vector{Float32}

    function VariantMap(
        chr::Vector{Int8},
        pos::Vector{UInt32},
        ref::Vector{Char},
        alt::Vector{Char},
        frq::Vector{Float32} = Float32[],
    )
        n = length(chr)
        length(pos) == n && length(ref) == n && length(alt) == n ||
            throw(ArgumentError("chr, pos, ref, and alt vectors must have the same length"))
        if !isempty(frq)
            length(frq) == n || throw(ArgumentError("frq must have the same length as other columns"))
        end
        new(chr, pos, ref, alt, frq)
    end
end

"""
    VariantMap(
        chr::Vector{<:Integer},
        pos::Vector{<:Integer},
        ref::Vector{Char},
        alt::Vector{Char};
        frq::Vector{<:Real} = Float32[],
    ) -> VariantMap

Construct a `VariantMap` from chromosome, position, and allele vectors.

Integer vectors for chromosome and position are converted to `Int8` and `UInt32` respectively.
If provided, allele frequency vector `frq` is converted to `Float32`.

# Arguments
- `chr::Vector{<:Integer}`: Chromosome numbers.
- `pos::Vector{<:Integer}`: Physical coordinates in base pairs.
- `ref::Vector{Char}`: Reference alleles.
- `alt::Vector{Char}`: Alternate alleles.

# Keywords
- `frq::Vector{<:Real}`: Optional alternate allele frequencies (default: `Float32[]`).

# Returns
- `VariantMap`: Initialized variant map.
"""
function VariantMap(
    chr::Vector{<:Integer},
    pos::Vector{<:Integer},
    ref::Vector{Char},
    alt::Vector{Char};
    frq::Vector{<:Real} = Float32[],
)
    VariantMap(
        Int8.(chr),
        UInt32.(pos),
        ref,
        alt,
        isempty(frq) ? Float32[] : Float32.(frq),
    )
end

Base.length(vm::VariantMap) = length(vm.pos)
Base.size(vm::VariantMap) = (length(vm.pos),)

@inline function Base.getindex(vm::VariantMap, i::Integer)
    @boundscheck checkbounds(vm.pos, i)
    return (
        chr = vm.chr[i],
        pos = vm.pos[i],
        ref = vm.ref[i],
        alt = vm.alt[i],
        frq = isempty(vm.frq) ? NaN32 : vm.frq[i],
    )
end

function Base.getindex(vm::VariantMap, idxs::AbstractVector{<:Integer})
    @boundscheck checkbounds(vm.pos, idxs)
    return VariantMap(
        vm.chr[idxs],
        vm.pos[idxs],
        vm.ref[idxs],
        vm.alt[idxs],
        isempty(vm.frq) ? Float32[] : vm.frq[idxs],
    )
end

function Base.show(io::IO, vm::VariantMap)
    print(io, "VariantMap with $(length(vm)) variants")
end

function Base.show(io::IO, ::MIME"text/plain", vm::VariantMap)
    println(io, "VariantMap with $(length(vm)) variants:")
    n = length(vm)
    if n <= 10
        for i = 1:n
            println(io, "  [$i] chr $(vm.chr[i]):$(vm.pos[i]) $(vm.ref[i])>$(vm.alt[i])" * (isempty(vm.frq) ? "" : " (frq: $(round(vm.frq[i], digits=4)))"))
        end
    else
        for i = 1:5
            println(io, "  [$i] chr $(vm.chr[i]):$(vm.pos[i]) $(vm.ref[i])>$(vm.alt[i])" * (isempty(vm.frq) ? "" : " (frq: $(round(vm.frq[i], digits=4)))"))
        end
        println(io, "  ⋮")
        for i = (n-4):n
            println(io, "  [$i] chr $(vm.chr[i]):$(vm.pos[i]) $(vm.ref[i])>$(vm.alt[i])" * (isempty(vm.frq) ? "" : " (frq: $(round(vm.frq[i], digits=4)))"))
        end
    end
end
