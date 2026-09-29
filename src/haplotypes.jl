"""
    struct Haplotype <: AbstractMatrix{Bool}

A compact, bit-packed matrix representation of phased haplotype data across genetic loci.

`Haplotype` organizes genetic data in a locus-major layout:
- Rows represent genetic loci (`1:nlc`).
- Columns represent phased haplotypes (`1:nhp`).
- For diploid individuals, each individual occupies two consecutive columns:
  individual `i` corresponds to haplotypes at column indices `2i - 1` and `2i`.

`Haplotype` implements the `AbstractMatrix{Bool}` interface (`size`, scalar indexing `h[locus, haplotype]`,
and assignment). Internally, data are stored in a `BitMatrix` whose row count is rounded up to a multiple
of 64 (`64 * cld(nlc, 64)`) to enable fast 64-bit word-level operations. Storage padding rows beyond
`nlc` are zero-masked and hidden from standard indexing.

# Fields
- `nlc::Int`: Number of genetic loci (logical rows; must be positive).
- `nhp::Int`: Number of haplotypes (columns; must be positive and even).
- `gt::BitMatrix`: Internal bit matrix of dimensions `(64 * cld(nlc, 64), nhp)`.

# Constructors
- `Haplotype(nlc::Int, nhp::Int, gt::BitMatrix)`: Low-level constructor wrapping an existing `BitMatrix`
  with matching dimensions; zero-masks padding bits beyond `nlc`.
- `Haplotype(nlc::Int, nhp::Int)`: Allocate an empty (all-false) `Haplotype` matrix.
- `Haplotype(gt::BitMatrix)`: Create a `Haplotype` by inferring `nlc = size(gt, 1)` and `nhp = size(gt, 2)`
  and copying the data into padded storage.

# Examples
```julia
using BnGStructs

# Create a haplotype matrix for 1,000 loci across 50 diploid individuals (100 haplotypes)
hps = Haplotype(1_000, 100)
size(hps) # (1000, 100)

# Set allele for locus 10 on the first haplotype of individual 1
hps[10, 1] = true
hps[10, 1] # true
```
"""
struct Haplotype <: AbstractMatrix{Bool}
    nlc::Int   # number of loci (rows logically = nlc, stored rows = 64 * cld(nlc, 64))
    nhp::Int   # number of haplotypes (columns in gt)
    gt::BitMatrix
    function Haplotype(nlc::Int, nhp::Int, gt::BitMatrix)
        nlc > 0 || error("nlc must be positive")
        (nhp > 0 && iseven(nhp)) || error("nhp must be positive and even")
        size(gt, 1) == 64 * cld(nlc, 64) || error("Number of rows of gt does not match nlc")
        size(gt, 2) == nhp || error("Number of columns of gt does not match nhp")
        # zero out padded rows beyond nlc
        _mask_last_word_columns!(gt, nlc)
        new(nlc, nhp, gt)
    end
end

"""
    Haplotype(nlc::Int, nhp::Int) -> Haplotype

Construct an all-zero (`false`) `Haplotype` object with `nlc` loci and `nhp` haplotypes.

# Arguments
- `nlc::Int`: Number of loci (must be positive).
- `nhp::Int`: Number of haplotypes (must be positive and even).

# Returns
- `Haplotype`: An initialized haplotype matrix backed by a zero-filled, 64-bit-aligned `BitMatrix`.
"""
function Haplotype(nlc::Int, nhp::Int)
    nlc > 0 || error("nlc must be positive")
    (nhp > 0 && iseven(nhp)) || error("nhp must be positive and even")
    gt = falses(64 * cld(nlc, 64), nhp) # as a chunk in BitMatrix is 64-bit UInt64
    return Haplotype(nlc, nhp, gt)
end

"""
    Haplotype(gt::BitMatrix) -> Haplotype

Construct a `Haplotype` from an existing `BitMatrix` `gt`.

The number of loci (`nlc`) and haplotypes (`nhp`) are inferred from `size(gt, 1)` and `size(gt, 2)`.
`nhp` must be positive and even. Only the valid rows (`1:nlc`) are copied into the internal 64-bit
chunk-aligned storage.

# Arguments
- `gt::BitMatrix`: Source bit matrix with dimensions `(nlc, nhp)`.

# Returns
- `Haplotype`: A new `Haplotype` instance containing the copied bits.
"""
function Haplotype(gt::BitMatrix)
    nlc, nhp = size(gt)
    th = Haplotype(nlc, nhp)
    th.gt[1:nlc, 1:nhp] = gt # copy only the valid part
    return th
end

function Base.show(io::IO, h::Haplotype)
    print(io, "Haplotype with $(h.nlc) loci and $(h.nhp) haplotypes")
end

# AbstractMatrix interface
Base.size(h::Haplotype) = (h.nlc, h.nhp)
Base.IndexStyle(::Type{<:Haplotype}) = IndexCartesian()

@inline function Base.getindex(h::Haplotype, i::Int, j::Int)
    @boundscheck checkbounds(h, i, j)
    @inbounds h.gt[i, j]
end

@inline function Base.setindex!(h::Haplotype, v, i::Int, j::Int)
    @boundscheck checkbounds(h, i, j)
    @inbounds h.gt[i, j] = Bool(v)
end
