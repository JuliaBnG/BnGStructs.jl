using BnGStructs
using Test
using Random

@testset "Haplotype to Genotype and back" begin
    nlc, nid = rand(500:1000), rand(500:1000)
    hps = Haplotype(nlc, 2 * nid)
    # fill in some random data
    rand!(hps.gt)
    hps = Haplotype(nlc, 2 * nid, hps.gt) # to mask out padded rows
    gt = hps.gt[1:nlc, 1:2:end] + hps.gt[1:nlc, 2:2:end]
    grm = gt'gt
    zz = gt * gt'

    id = hap2id(hps)
    gt = id.gt[1:nid, 1:2:end] + id.gt[1:nid, 2:2:end]
    grm2 = gt * gt'
    zz2 = gt'gt

    hps2 = id2hap(id)
    gt = hps2.gt[1:nlc, 1:2:end] + hps2.gt[1:nlc, 2:2:end]
    grm3 = gt'gt
    zz3 = gt * gt'

    @test grm == grm2
    @test zz == zz2
    @test grm == grm3
    @test zz == zz3
    @test hps.gt == hps2.gt
end