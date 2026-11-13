"""
    struct VariantMap

Coordinates and alleles of the variants of a genome, as parallel vectors.

# Coordinate convention

`pos` is the **genome-wide** base-pair coordinate: chromosomes are laid end
to end in the order of the species, so a variant at base `b` of chromosome
`c` has `pos = cbp(sp)[c-1] + b`. This is the coordinate FisherWright
simulates on and the one meiosis in BnGsim searches. The inner constructor
enforces what follows from it:

- `pos` is strictly increasing (sorted, no two variants at one site);
- `chr` is non-decreasing.

Use [`chrpos`](@ref) for the within-chromosome coordinate,
`VariantMap(chr, chrpos, sp)` to build a map from within-chromosome
coordinates, and [`check_coordinates`](@ref) to verify a map against a
species.

# Fields
- `chr::Vector{Int8}`: Chromosome numbers, 1-based.
- `pos::Vector{UInt32}`: Genome-wide base-pair coordinates (see above).
- `ref::Vector{Char}`: Reference alleles (haplotype bit 0).
- `alt::Vector{Char}`: Alternate alleles (haplotype bit 1).
- `frq::Vector{Float32}`: Frequency of the alternate allele; empty if not recorded.
- `morgan::Vector{Float64}`: Genetic position in Morgans **within the
  chromosome**, non-decreasing within each chromosome; empty if the map is
  the uniform one implied by `M(sp)`.

# Constructors
- `VariantMap(chr, pos, ref, alt[, frq]; frq, morgan)`: from genome-wide coordinates.
- `VariantMap(pos, sp::Species; ref, alt, frq, morgan)`: `chr` derived from `pos`.
- `VariantMap(chr, chrpos, sp::Species; ref, alt, frq, morgan)`: from
  within-chromosome coordinates.

# Indexing
- `length(vm)`, `size(vm)`.
- `vm[i]`: `(chr, pos, ref, alt, frq)` of variant `i`; `frq` is `NaN32` if not recorded.
- `vm[idxs]`: the sub-map at `idxs`, which must be increasing (a range, a
  sorted index vector or a `LocusSet`'s `loci`).

# Examples
```julia
using BnGStructs

sp = GenericSpecies("Toy", 10, [5_000, 5_000])
# Variants at bases 1000 and 2000 of chromosome 1 and base 1500 of chromosome 2
vm = VariantMap([1, 1, 2], [1000, 2000, 1500], sp)
vm.pos           # UInt32[1000, 2000, 6500]
chrpos(vm, sp)   # UInt32[1000, 2000, 1500]
sub = vm[[1, 3]] # VariantMap with 2 variants
```
"""
struct VariantMap
    chr::Vector{Int8}
    pos::Vector{UInt32}
    ref::Vector{Char}
    alt::Vector{Char}
    frq::Vector{Float32}
    morgan::Vector{Float64}

    function VariantMap(
        chr::Vector{Int8},
        pos::Vector{UInt32},
        ref::Vector{Char},
        alt::Vector{Char},
        frq::Vector{Float32},
        morgan::Vector{Float64},
    )
        n = length(chr)
        length(pos) == n && length(ref) == n && length(alt) == n ||
            throw(ArgumentError("chr, pos, ref, and alt vectors must have the same length"))
        isempty(frq) || length(frq) == n ||
            throw(ArgumentError("frq must have the same length as other columns"))
        isempty(morgan) || length(morgan) == n ||
            throw(ArgumentError("morgan must have the same length as other columns"))
        @inbounds for i = 2:n
            pos[i-1] < pos[i] || throw(ArgumentError(
                "pos must be strictly increasing genome-wide coordinates; " *
                "pos[$(i - 1)] = $(Int(pos[i-1])) and pos[$i] = $(Int(pos[i])). " *
                "For within-chromosome positions use VariantMap(chr, chrpos, species).",
            ))
            chr[i-1] <= chr[i] ||
                throw(ArgumentError("chr must be non-decreasing; chr[$(i - 1)] = $(chr[i-1]) and chr[$i] = $(chr[i])"))
            isempty(morgan) || chr[i-1] != chr[i] || morgan[i-1] <= morgan[i] ||
                throw(ArgumentError("morgan must be non-decreasing within a chromosome (variant $i)"))
        end
        all(>(0), chr) || throw(ArgumentError("chromosome numbers are 1-based"))
        isempty(morgan) || all(m -> isfinite(m) && m >= 0, morgan) ||
            throw(ArgumentError("morgan must be finite and non-negative"))
        new(chr, pos, ref, alt, frq, morgan)
    end
end

function VariantMap(
    chr::AbstractVector{<:Integer},
    pos::AbstractVector{<:Integer},
    ref::AbstractVector{Char},
    alt::AbstractVector{Char},
    frq_positional::AbstractVector{<:Real} = Float32[];
    frq::AbstractVector{<:Real} = frq_positional,
    morgan::AbstractVector{<:Real} = Float64[],
)
    VariantMap(
        convert(Vector{Int8}, chr),
        convert(Vector{UInt32}, pos),
        convert(Vector{Char}, ref),
        convert(Vector{Char}, alt),
        convert(Vector{Float32}, frq),
        convert(Vector{Float64}, morgan),
    )
end

"""
    VariantMap(pos, sp::Species; ref, alt, frq = Float32[], morgan = Float64[])

Build a map from genome-wide coordinates `pos`, taking the chromosome of
each variant from the chromosome ends of `sp`. `ref` and `alt` default to
the placeholders `'A'` and `'C'`.
"""
function VariantMap(
    pos::AbstractVector{<:Integer},
    sp::Species;
    ref::AbstractVector{Char} = fill('A', length(pos)),
    alt::AbstractVector{Char} = fill('C', length(pos)),
    frq::AbstractVector{<:Real} = Float32[],
    morgan::AbstractVector{<:Real} = Float64[],
)
    ends = cbp(sp)
    isempty(pos) || 1 <= minimum(pos) && maximum(pos) <= ends[end] ||
        throw(ArgumentError("pos must lie in 1:$(Int(ends[end])), the genome of $(sp.name)"))
    chr = [Int8(searchsortedfirst(ends, p)) for p in pos]
    VariantMap(chr, pos, ref, alt; frq = frq, morgan = morgan)
end

"""
    VariantMap(chr, chrpos, sp::Species; ref, alt, frq = Float32[], morgan = Float64[])

Build a map from chromosome numbers and **within-chromosome** base-pair
positions, converting to the genome-wide coordinate with the chromosome
lengths of `sp`. Variants must be ordered by chromosome and position.
"""
function VariantMap(
    chr::AbstractVector{<:Integer},
    chrpos::AbstractVector{<:Integer},
    sp::Species;
    ref::AbstractVector{Char} = fill('A', length(chr)),
    alt::AbstractVector{Char} = fill('C', length(chr)),
    frq::AbstractVector{<:Real} = Float32[],
    morgan::AbstractVector{<:Real} = Float64[],
)
    length(chr) == length(chrpos) ||
        throw(ArgumentError("chr and chrpos must have the same length"))
    ends = cbp(sp)
    len = sp.chromosome
    pos = Vector{UInt32}(undef, length(chr))
    for i in eachindex(chr)
        c, b = chr[i], chrpos[i]
        1 <= c <= length(len) ||
            throw(ArgumentError("variant $i is on chromosome $c; $(sp.name) has $(length(len))"))
        1 <= b <= len[c] ||
            throw(ArgumentError("variant $i is at base $b of chromosome $c, which has $(Int(len[c])) bp"))
        pos[i] = (c == 1 ? UInt32(0) : ends[c-1]) + UInt32(b)
    end
    VariantMap(chr, pos, ref, alt; frq = frq, morgan = morgan)
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
        isempty(vm.morgan) ? Float64[] : vm.morgan[idxs],
    )
end

"""
    chrpos(vm::VariantMap, sp::Species) -> Vector{UInt32}
    chrpos(vm::VariantMap, sp::Species, i::Integer) -> UInt32

Within-chromosome base-pair position of every variant, or of variant `i`:
`vm.pos` minus the end of the preceding chromosome of `sp`.
"""
function chrpos(vm::VariantMap, sp::Species, i::Integer)
    c = vm.chr[i]
    return c == 1 ? vm.pos[i] : vm.pos[i] - cbp(sp)[c-1]
end

function chrpos(vm::VariantMap, sp::Species)
    ends = cbp(sp)
    return [c == 1 ? p : p - ends[c-1] for (c, p) in zip(vm.chr, vm.pos)]
end

"""
    check_coordinates(vm::VariantMap, sp::Species) -> true

Verify that `vm` is on the genome of `sp`: every `pos` lies within the
genome, and `chr` is the chromosome that contains it. Throws an
`ArgumentError` naming the first variant that fails. A map built from
within-chromosome positions by mistake fails here unless it has a single
chromosome.
"""
function check_coordinates(vm::VariantMap, sp::Species)
    ends = cbp(sp)
    for i in eachindex(vm.pos)
        p = vm.pos[i]
        1 <= p <= ends[end] || throw(ArgumentError(
            "variant $i: pos = $(Int(p)) is outside the genome of $(sp.name) (1:$(Int(ends[end])))"))
        c = searchsortedfirst(ends, p)
        c == vm.chr[i] || throw(ArgumentError(
            "variant $i: pos = $(Int(p)) lies on chromosome $c of $(sp.name), but chr = $(vm.chr[i]). " *
            "pos must be the genome-wide coordinate."))
    end
    return true
end

function Base.show(io::IO, vm::VariantMap)
    print(io, "VariantMap with $(length(vm)) variants")
end

function _show_variant(io::IO, vm::VariantMap, i::Int)
    print(io, "  [$i] chr $(vm.chr[i]) @ $(Int(vm.pos[i])) $(vm.ref[i])>$(vm.alt[i])")
    isempty(vm.frq) || print(io, " (frq: $(round(vm.frq[i], digits=4)))")
    isempty(vm.morgan) || print(io, " (M: $(round(vm.morgan[i], digits=4)))")
    println(io)
end

function Base.show(io::IO, ::MIME"text/plain", vm::VariantMap)
    println(io, "VariantMap with $(length(vm)) variants (genome-wide bp):")
    n = length(vm)
    if n <= 10
        foreach(i -> _show_variant(io, vm, i), 1:n)
    else
        foreach(i -> _show_variant(io, vm, i), 1:5)
        println(io, "  ⋮")
        foreach(i -> _show_variant(io, vm, i), (n-4):n)
    end
end
