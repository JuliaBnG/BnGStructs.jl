"""
    hap2id(hps::Haplotype)
This is to transform a `Haplotype` object to a `Genotype` object.

Given the SNP haplotypes `hps` of type `Haplotype`, this function transposes its
`gt` field (of size `nlc × 2nid`) to a `BitMatrix` of size `n × 2nlc`, where `n`
is the minimum number that is a multiple of 64 and is no less than `nid`. The
resulting `BitMatrix` is used to construct and return a `Genotype` object.

Note this function deems the computer is of 64-bit architecture.
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
    id2hap(g::Genotype)
This is to transform a `Genotype` object back to a `Haplotype` object,
performing the reverse operation of `hap2id`.

Given a `Genotype` object `g`, this function transposes its `gt` field (of size
`nid × 2nlc`) to a `BitMatrix` of size `m × 2nid`, where `m` is the minimum
number that is a multiple of 64 and is no less than `nlc`. The resulting
`BitMatrix` is used to construct and return a `Haplotype` object.

Note this function deems the computer is of 64-bit architecture.
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
