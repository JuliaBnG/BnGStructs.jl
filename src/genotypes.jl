"""
    struct Genotype <: AbstractMatrix{Bool}

A compact, bit-packed matrix representation of phased individual genotype data.

`Genotype` organizes genetic data in an individual-major layout:
- Rows represent individuals (`1:nid`).
- Columns represent phased alleles across all loci (`1:nas`).
- For diploid individuals with `nlc` loci, `nas = 2 * nlc`. For locus `l`, the two homologous
  alleles reside at column indices `2l - 1` and `2l`.

`Genotype` implements the `AbstractMatrix{Bool}` interface (`size`, scalar indexing `g[individual, allele]`,
and assignment). Internally, data are stored in a `BitMatrix` whose row count is rounded up to a multiple
of 64 (`64 * cld(nid, 64)`) to enable fast 64-bit word-level operations along individuals. Storage padding rows
beyond `nid` are zero-masked and hidden from standard indexing.

# Fields
- `nid::Int`: Number of individuals (logical rows; must be positive).
- `nas::Int`: Number of alleles (columns; must be positive and even, equal to `2 * n_loci`).
- `gt::BitMatrix`: Internal bit matrix of dimensions `(64 * cld(nid, 64), nas)`.

# Constructors
- `Genotype(nid::Int, nas::Int, gt::BitMatrix)`: Low-level constructor wrapping an existing `BitMatrix`
  with matching dimensions; zero-masks padding bits beyond `nid`.
- `Genotype(nid::Int, nas::Int)`: Allocate an empty (all-false) `Genotype` matrix.
- `Genotype(gt::BitMatrix)`: Create a `Genotype` by inferring `nid = size(gt, 1)` and `nas = size(gt, 2)`
  and copying the data into padded storage.

# Examples
```julia
using BnGStructs

# Create a genotype matrix for 50 individuals across 1,000 diploid loci (2,000 alleles)
gt = Genotype(50, 2_000)
size(gt) # (50, 2000)

# Set the first allele of individual 1 at locus 5 (allele index 2*5 - 1 = 9)
gt[1, 9] = true
gt[1, 9] # true
```
"""
struct Genotype <: AbstractMatrix{Bool}
    nid::Int   # number of individuals
    nas::Int   # number of alleles
    gt::BitMatrix
    function Genotype(nid::Int, nas::Int, gt::BitMatrix)
        nid > 0 || error("nid must be positive")
        (nas > 0 && iseven(nas)) || error("nas must be positive and even")
        size(gt, 1) == 64 * cld(nid, 64) || error("Number of rows of gt does not match nid")
        size(gt, 2) == nas || error("Number of columns of gt does not match nas")
        # zero out padded rows beyond nid
        _mask_last_word_columns!(gt, nid)
        new(nid, nas, gt)
    end
end

"""
    Genotype(nid::Int, nas::Int) -> Genotype

Construct an all-zero (`false`) `Genotype` object with `nid` individuals and `nas` alleles.

# Arguments
- `nid::Int`: Number of individuals (must be positive).
- `nas::Int`: Number of alleles (must be positive and even, typically `2 * nlc`).

# Returns
- `Genotype`: An initialized genotype matrix backed by a zero-filled, 64-bit-aligned `BitMatrix`.
"""
function Genotype(nid::Int, nas::Int)
    (nas > 0 && iseven(nas)) || error("nas must be positive even")
    nid > 0 || error("nid must be positive")
    gt = falses(64 * cld(nid, 64), nas) # as a chunk in BitMatrix is 64-bit UInt64
    return Genotype(nid, nas, gt)
end

"""
    Genotype(gt::BitMatrix) -> Genotype

Construct a `Genotype` from an existing `BitMatrix` `gt`.

The number of individuals (`nid`) and alleles (`nas`) are inferred from `size(gt, 1)` and `size(gt, 2)`.
`nas` must be positive and even. Only the valid rows (`1:nid`) are copied into the internal 64-bit
chunk-aligned storage.

# Arguments
- `gt::BitMatrix`: Source bit matrix with dimensions `(nid, nas)`.

# Returns
- `Genotype`: A new `Genotype` instance containing the copied bits.
"""
function Genotype(gt::BitMatrix)
    nid, nas = size(gt)
    tg = Genotype(nid, nas)
    tg.gt[1:nid, 1:nas] = gt[1:nid, 1:nas] # copy only the valid part
    return tg
end

function Base.show(io::IO, g::Genotype)
    print(io, "Genotype with $(g.nid) individuals and $(g.nas) alleles")
end

# AbstractMatrix interface
Base.size(g::Genotype) = (g.nid, g.nas)
Base.IndexStyle(::Type{<:Genotype}) = IndexCartesian()

@inline function Base.getindex(g::Genotype, i::Int, j::Int)
    @boundscheck checkbounds(g, i, j)
    @inbounds g.gt[i, j]
end

@inline function Base.setindex!(g::Genotype, v, i::Int, j::Int)
    @boundscheck checkbounds(g, i, j)
    @inbounds g.gt[i, j] = Bool(v)
end
