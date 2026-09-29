"""
    hap2id(hps::Haplotype) -> Genotype

Convert a locus-major `Haplotype` matrix into an individual-major `Genotype` matrix.

This operation transposes the underlying genetic bits:
- Input `hps` has dimensions `(nlc, 2 * nid)`.
- Output `Genotype` has dimensions `(nid, 2 * nlc)`.

For each individual `i ∈ 1:nid` and locus `l ∈ 1:nlc`:
- `hps[l, 2i - 1]` maps to `gt[i, 2l - 1]`.
- `hps[l, 2i]` maps to `gt[i, 2l]`.

The transposition operates directly on 64-bit word chunks in parallel using
multi-threading (`Threads.@threads`).

# Arguments
- `hps::Haplotype`: Source haplotype matrix.

# Returns
- `Genotype`: Transposed genotype matrix.

# Examples
```julia
using BnGStructs

hps = Haplotype(500, 200) # 500 loci, 100 individuals
gt = hap2id(hps)
size(gt) # (100, 1000)
```

See also: [`id2hap`](@ref).
"""
function hap2id(hps::Haplotype)
    nlc, nhp = hps.nlc, hps.nhp
    nid = nhp ÷ 2
    nas = 2nlc
    
    # The Genotype constructor creates a correctly padded BitMatrix of falses
    gt = Genotype(nid, nas)

    _transpose_gt!(gt.gt.chunks, hps.gt, nlc, nid)
    
    return gt
end

"""
    id2hap(g::Genotype) -> Haplotype

Convert an individual-major `Genotype` matrix back into a locus-major `Haplotype` matrix.

Performs the inverse transformation of [`hap2id`](@ref):
- Input `g` has dimensions `(nid, 2 * nlc)`.
- Output `Haplotype` has dimensions `(nlc, 2 * nid)`.

For each individual `i ∈ 1:nid` and locus `l ∈ 1:nlc`:
- `g[i, 2l - 1]` maps to `hps[l, 2i - 1]`.
- `g[i, 2l]` maps to `hps[l, 2i]`.

The transposition operates directly on 64-bit word chunks in parallel using
multi-threading (`Threads.@threads`).

# Arguments
- `g::Genotype`: Source genotype matrix.

# Returns
- `Haplotype`: Transposed haplotype matrix.

# Examples
```julia
using BnGStructs

gt = Genotype(100, 1_000) # 100 individuals, 500 loci
hps = id2hap(gt)
size(hps) # (500, 200)
```

See also: [`hap2id`](@ref).
"""
function id2hap(g::Genotype)
    nid, nas = g.nid, g.nas
    nlc = nas ÷ 2
    nhp = 2nid

    # The Haplotype constructor creates a correctly padded BitMatrix of falses
    hps = Haplotype(nlc, nhp)

    _transpose_gt!(hps.gt.chunks, g.gt, nid, nlc)

    return hps
end
