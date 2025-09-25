"""
    _mask_last_word_columns!(bm::BitMatrix, nrow::Int)
Zero out the bits in the last 64-bit word of each column of the BitMatrix `bm`
beyond the `nrow`-th row. This is useful when the number of logical rows is not a
multiple of 64, and the BitMatrix has been allocated with extra rows to accommodate
the storage in 64-bit words.
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
    _transpose_gt!(dest_chunks, src_gt, outer_dim, inner_dim)

A private helper function to transpose genotype data between haplotype-major and
individual-major layouts. It operates directly on the `chunks` of the underlying
`BitMatrix` for performance.

# Arguments
- `dest_chunks`: The `chunks` vector of the destination `BitMatrix`.
- `src_gt`: The source `BitMatrix`.
- `outer_dim`: The dimension to parallelize over (e.g., `nlc` or `nid`).
- `inner_dim`: The other dimension.
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
