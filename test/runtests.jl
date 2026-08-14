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

    # Test AbstractMatrix indexing on Haplotype and Genotype
    @test size(hps) == (nlc, 2 * nid)
    @test hps[1, 1] isa Bool
    @test hps[1, 1] == hps.gt[1, 1]
    hps[1, 1] = !hps[1, 1]
    @test hps[1, 1] == hps.gt[1, 1]

    @test size(id) == (nid, 2 * nlc)
    @test id[1, 1] isa Bool
    @test id[1, 1] == id.gt[1, 1]
    id[1, 1] = !id[1, 1]
    @test id[1, 1] == id.gt[1, 1]
end

@testset "Species constructors and accessors" begin
    cattle = Cattle(1000; M=80_000_000)
    @test cattle isa Species
    @test cattle.name == "BosTau"
    @test cattle.nid == 1000
    @test name(cattle) == "BosTau"
    @test nid(cattle) == 1000
    @test chromosome(cattle) == cattle.chromosome
    @test M(cattle) == 80_000_000
    @test length(cattle) == 29
    @test total_bp(cattle) == sum(cattle.chromosome)
    @test cbp(cattle) == cumsum(cattle.chromosome)
    @test cattle.M == 80_000_000

    sheep = Sheep(500)
    @test sheep isa Sheep
    @test sheep.M == 100_000_000

    chicken = Chicken(200; M=90_000_000)
    @test chicken isa Chicken
    @test chicken.M == 90_000_000

    pig = Pig(300)
    @test pig isa Pig
    goat = Goat(150)
    @test goat isa Goat
    horse = Horse(100)
    @test horse isa Horse
    rabbit = Rabbit(400)
    @test rabbit isa Rabbit
    cat = Cat(250)
    @test cat isa Cat
    dog = Dog(350)
    @test dog isa Dog

    gen = GenericSpecies("Custom", 500, [10_000_000, 20_000_000]; M=50_000_000)
    @test gen isa GenericSpecies
    @test gen.name == "Custom"
    @test gen.nid == 500
    @test length(gen) == 2
    @test total_bp(gen) == 30_000_000
    @test cbp(gen) == UInt32[10_000_000, 30_000_000]
    @test gen.M == 50_000_000

    @test_throws ErrorException GenericSpecies("Bad", -1, [1000])
    @test_throws ErrorException GenericSpecies("Bad", 10, Int[])
    @test_throws ErrorException GenericSpecies("Bad", 10, [0])
    @test_throws ErrorException GenericSpecies("Bad", 10, [1000]; M=0)

    for constructor in (Cattle, Sheep, Chicken, Pig, Goat, Horse, Rabbit, Cat, Dog)
        @test_throws ErrorException constructor(-1)
        @test_throws ErrorException constructor(100; M=0)
    end
end

@testset "Trait definitions" begin
    using Distributions: Normal, Gamma
    t1 = Trait("MilkYield"; h²=0.3, σₐ=2.5, da=Normal(0, 2.5))
    @test t1 isa BnGStructs.aTrait
    @test t1.name == "MilkYield"
    @test t1.h² == 0.3
    @test t1.da isa Normal

    t2 = Trait("DiseaseRisk", [0.7, 0.2, 0.1]; h²=0.2)
    @test t2 isa BnGStructs.tTrait
    @test t2.name == "DiseaseRisk"
    @test length(t2.threshold) == 2

    weights = [0.7, 0.2, 0.1]
    original_weights = copy(weights)
    @test Trait("IntegerWeights", [1, 2]) isa BnGStructs.tTrait
    Trait("UnchangedWeights", weights)
    @test weights == original_weights
    @test_throws ErrorException Trait("OneCategory", [1.0])
    @test_throws ErrorException Trait("InfiniteWeight", [Inf, 1.0])
    @test_throws ErrorException Trait("InfiniteAge"; age=Inf)
    @test_throws ErrorException Trait("InfiniteMean"; μ=Inf)
end

@testset "Container display" begin
    @test sprint(show, Haplotype(10, 2)) == "Haplotype with 10 loci and 2 haplotypes"
    @test sprint(show, Genotype(10, 2)) == "Genotype with 10 individuals and 2 alleles"
end

@testset "VariantMap and LocusSet" begin
    chrs = Int8[1, 1, 1, 2, 2]
    pos = UInt32[100, 200, 300, 150, 250]
    ref = ['A', 'C', 'G', 'T', 'A']
    alt = ['G', 'T', 'A', 'C', 'G']
    frq = Float32[0.1, 0.2, 0.5, 0.4, 0.3]

    vm = VariantMap(chrs, pos, ref, alt, frq)
    @test length(vm) == 5
    @test size(vm) == (5,)
    @test vm[1].pos == 100
    @test vm[1].ref == 'A'
    @test vm[1].frq == 0.1f0

    sub_vm = vm[[1, 3, 5]]
    @test length(sub_vm) == 3
    @test sub_vm.pos == UInt32[100, 300, 250]

    # Locus sets (Chip markers, QTLs, Reference loci)
    chip = LocusSet("50kChip", [1, 2, 4])
    qtl_set = LocusSet("10kQTL", [2, 3, 5])
    ref_set = LocusSet("RefPanel", [1, 5])

    @test length(chip) == 3
    @test 2 in chip
    @test !(3 in chip)
    @test length(qtl_set) == 3
    @test 5 in qtl_set
    @test sprint(show, chip) == "LocusSet \"50kChip\" with 3 loci"

    @test_throws ArgumentError VariantMap([1, 2], [100], ['A', 'C'], ['G', 'T'])
    @test_throws ArgumentError LocusSet("Bad", [0, 1])
    @test_throws ArgumentError LocusSet("Unsorted", [2, 1])
    @test_throws ArgumentError LocusSet("Duplicate", [1, 1])
end

@testset "TraitQTL, MultiTraitQTL and TBV" begin
    # 4 loci, 2 individuals (4 haplotypes)
    hps = Haplotype(4, 4)
    # Individual 1: hap1 = [1, 0, 1, 0], hap2 = [1, 1, 0, 0] -> dosages = [2, 1, 1, 0]
    # Individual 2: hap3 = [0, 0, 0, 0], hap4 = [0, 1, 1, 1] -> dosages = [0, 1, 1, 1]
    hps.gt[1, 1] = true; hps.gt[3, 1] = true
    hps.gt[1, 2] = true; hps.gt[2, 2] = true
    hps.gt[2, 4] = true; hps.gt[3, 4] = true; hps.gt[4, 4] = true

    # Trait 1: QTL at loci 1 and 2 with effects [0.5, -1.0]
    qtl1 = TraitQTL("Trait1", [1, 2], [0.5, -1.0])
    @test length(qtl1) == 2
    tbv1 = tbv(hps, qtl1)
    # Ind 1: 2*0.5 + 1*(-1.0) = 0.0
    # Ind 2: 0*0.5 + 1*(-1.0) = -1.0
    @test tbv1 ≈ [0.0, -1.0]

    # Same with Genotype
    gt = hap2id(hps)
    @test tbv(gt, qtl1) ≈ [0.0, -1.0]

    # Trait 2: QTL at loci 2 and 3 with additive [2.0, 3.0] and dominance on locus 2 [0.5, 0.0]
    qtl2 = TraitQTL("Trait2", [2, 3], [2.0, 3.0]; dominance=[0.5, 0.0])
    tbv2 = tbv(hps, qtl2)
    # Ind 1: locus2 (dosage 1) -> 1*2.0 + 0.5 = 2.5; locus3 (dosage 1) -> 1*3.0 = 3.0. Total = 5.5
    # Ind 2: locus2 (dosage 1) -> 1*2.0 + 0.5 = 2.5; locus3 (dosage 1) -> 1*3.0 = 3.0. Total = 5.5
    @test tbv2 ≈ [5.5, 5.5]
    @test tbv(gt, qtl2) ≈ [5.5, 5.5]

    # MultiTraitQTL combining Trait 1 and Trait 2
    mqtl = MultiTraitQTL([qtl1, qtl2])
    @test length(mqtl.trait_names) == 2
    @test mqtl.loci == Int32[1, 2, 3] # union of QTL loci
    @test size(mqtl.additive) == (3, 2)

    tbv_multi = tbv(hps, mqtl)
    @test size(tbv_multi) == (2, 2)
    @test tbv_multi[:, 1] ≈ tbv1
    @test tbv_multi[:, 2] ≈ tbv2

    @test tbv(gt, mqtl) ≈ tbv_multi
    @test_throws ArgumentError TraitQTL("Duplicate", [1, 1], [0.1, 0.2])
end