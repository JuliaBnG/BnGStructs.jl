"""
    _mask_last_word_columns!(bm::BitMatrix, nrow::Int)

Zero out bits beyond the `nrow`-th logical row in the final 64-bit word (`UInt64` chunk)
of each column of `bm`.

In Julia's `BitMatrix`, columns are stored as sequences of 64-bit words (`UInt64`).
When the number of logical rows `nrow` is not an exact multiple of 64, extra bits in
the final word of each column serve as storage padding. This in-place helper clears
those excess padding bits to ensure consistent bitwise comparisons and hashing.

# Arguments
- `bm::BitMatrix`: The bit matrix whose terminal chunks will be masked.
- `nrow::Int`: The number of logical rows.
"""
function _mask_last_word_columns!(bm::BitMatrix, nrow::Int)
    m, n = size(bm)
    if nrow < m
        # Keep only the valid (lower) bits in the final 64-bit word of each column.
        # Number of valid bits in the last word:
        k = nrow & 0x3f         # same as nrow % 64, but faster
        @assert k != 0          # because nrow < m implies there is padding
        mask = (UInt64(1) << k) - UInt64(1)
        words_per_col = m >>> 6 # divide by 64
        chunks = bm.chunks
        @inbounds for col = 0:(n-1)
            idx = col * words_per_col + words_per_col
            chunks[idx] &= mask
        end
    end
end

"""
    _transpose_gt!(dest_chunks::Vector{UInt64}, src_gt::BitMatrix, outer_dim::Int, inner_dim::Int)

Transpose diploid genotype bits between locus-major (`Haplotype`) and individual-major (`Genotype`) layouts.

This internal function writes directly into the preallocated 64-bit `UInt64` chunks of the destination
`BitMatrix`. Multi-threading is parallelized across `outer_dim` (`Threads.@threads`) where each outer
index corresponds to two adjacent columns in the destination matrix (representing the two homologous alleles
or haplotypes per individual).

# Arguments
- `dest_chunks::Vector{UInt64}`: The chunk vector (`bm.chunks`) of the destination `BitMatrix`.
- `src_gt::BitMatrix`: The source genotype/haplotype `BitMatrix`.
- `outer_dim::Int`: The outer dimension size to parallelize over (`nlc` when transposing haplotype to genotype; `nid` when transposing genotype to haplotype).
- `inner_dim::Int`: The inner dimension size (`nid` when transposing haplotype to genotype; `nlc` when transposing genotype to haplotype).
"""
function _transpose_gt!(
    dest_chunks::Vector{UInt64},
    src_gt::BitMatrix,
    outer_dim::Int,
    inner_dim::Int
)
    chunks_per_dest_col = cld(inner_dim, 64)

    @inbounds Threads.@threads for outer_idx = 1:outer_dim
        dest_col1_start_chunk = (2outer_idx - 2) * chunks_per_dest_col
        dest_col2_start_chunk = (2outer_idx - 1) * chunks_per_dest_col

        for inner_chunk_idx = 1:chunks_per_dest_col
            inner_start = (inner_chunk_idx - 1) * 64 + 1
            
            chunk1 = UInt64(0)
            chunk2 = UInt64(0)

            for k = 0:63
                inner_idx = inner_start + k
                if inner_idx <= inner_dim
                    h1 = src_gt[outer_idx, 2 * inner_idx - 1]
                    h2 = src_gt[outer_idx, 2 * inner_idx]
                    
                    chunk1 |= UInt64(h1) << k
                    chunk2 |= UInt64(h2) << k
                end
            end
            
            dest_chunks[dest_col1_start_chunk + inner_chunk_idx] = chunk1
            dest_chunks[dest_col2_start_chunk + inner_chunk_idx] = chunk2
        end
    end
end
