const _BNGH1_MAGIC = UInt8['B', 'N', 'G', 'H', '1', 0, 0, 0]

"""
    write_haplotypes(path, hap::Haplotype, pos) -> path
    write_haplotypes(path, hap::Haplotype, vm::VariantMap) -> path

Write phased haplotypes and their genome-wide positions in the `BNGH1`
format. `pos` has one coordinate per locus of `hap`; with a
[`VariantMap`](@ref), `vm.pos` is written.

# Format (little-endian)

| Bytes | Content |
| :--- | :--- |
| 8 | magic `"BNGH1\\0\\0\\0"` |
| 8 + 8 | `Int64` number of loci, `Int64` number of haplotypes |
| 4 × nloci | `UInt32` genome-wide position of each locus |
| nloci × ⌈nhap / 8⌉ | one record per locus: bit `h - 1` (least significant first) is the allele of haplotype `h` |

Chromosome, alleles and frequencies are not stored; rebuild the map with
`VariantMap(pos, species)`. The format is fixed: a change of layout gets a
new magic.
"""
function write_haplotypes(path::AbstractString, hap::Haplotype, pos::AbstractVector{<:Integer})
    nlc, nhp = hap.nlc, hap.nhp
    length(pos) == nlc || throw(DimensionMismatch("pos needs one coordinate per locus ($nlc)"))
    nb = cld(nhp, 8)
    rec = Matrix{UInt8}(undef, nb, nlc)
    fill!(rec, 0x00)
    gt = hap.gt
    @inbounds for h = 1:nhp
        b, s = ((h - 1) >>> 3) + 1, (h - 1) & 7
        for l = 1:nlc
            gt[l, h] && (rec[b, l] |= 0x01 << s)
        end
    end
    open(path, "w") do io
        write(io, _BNGH1_MAGIC)
        write(io, htol(Int64(nlc)), htol(Int64(nhp)))
        write(io, htol.(convert(Vector{UInt32}, pos)))
        write(io, rec)
    end
    return path
end
write_haplotypes(path::AbstractString, hap::Haplotype, vm::VariantMap) =
    write_haplotypes(path, hap, vm.pos)

"""
    read_haplotypes(path) -> (hap::Haplotype, pos::Vector{UInt32})

Read a `BNGH1` file written by [`write_haplotypes`](@ref) or by another
program (see there for the layout). Loci come back in file order with
their positions as stored; nothing is sorted or filtered, so a file from
another simulator may hold monomorphic loci or unsorted positions.
"""
function read_haplotypes(path::AbstractString)
    open(path) do io
        magic = read(io, 8)
        magic == _BNGH1_MAGIC || throw(ArgumentError("$path is not a BNGH1 file"))
        nlc, nhp = Int(ltoh(read(io, Int64))), Int(ltoh(read(io, Int64)))
        nlc > 0 || throw(ArgumentError("$path: number of loci $nlc must be positive"))
        (nhp > 0 && iseven(nhp)) ||
            throw(ArgumentError("$path: number of haplotypes $nhp must be positive and even"))
        pos = Vector{UInt32}(undef, nlc)
        read!(io, pos)
        pos .= ltoh.(pos)
        nb = cld(nhp, 8)
        rec = Matrix{UInt8}(undef, nb, nlc)
        read!(io, rec)
        eof(io) || throw(ArgumentError("$path has trailing bytes after the haplotype records"))
        hap = Haplotype(nlc, nhp)
        gt = hap.gt
        Threads.@threads for h = 1:nhp
            b, s = ((h - 1) >>> 3) + 1, (h - 1) & 7
            @inbounds for l = 1:nlc
                gt[l, h] = (rec[b, l] >>> s) & 0x01 == 0x01
            end
        end
        return hap, pos
    end
end
