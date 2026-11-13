using BnGStructs
using Test
using Random
using Aqua

@testset "Aqua" begin
    Aqua.test_all(BnGStructs)
end

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
    @test BnGStructs.name(cattle) == "BosTau"
    @test !isdefined(@__MODULE__, :name) && !isdefined(@__MODULE__, :M)   # no longer exported
    @test nid(cattle) == 1000
    @test chromosome(cattle) == cattle.chromosome
    @test BnGStructs.M(cattle) == 80_000_000
    @test length(cattle) == 29
    @test total_bp(cattle) == sum(cattle.chromosome)
    @test cbp(cattle) == cumsum(cattle.chromosome)
    @test cattle.M == 80_000_000

    sheep = Sheep(500)
    @test sheep isa Species
    @test sheep.M == 100_000_000

    chicken = Chicken(200; M=90_000_000)
    @test chicken isa Species
    @test chicken.M == 90_000_000

    pig = Pig(300)
    @test pig isa Species
    goat = Goat(150)
    @test goat isa Species
    horse = Horse(100)
    @test horse isa Species
    rabbit = Rabbit(400)
    @test rabbit isa Species
    cat = Cat(250)
    @test cat isa Species
    dog = Dog(350)
    @test dog isa Species

    gen = GenericSpecies("Custom", 500, [10_000_000, 20_000_000]; M=50_000_000)
    @test gen isa Species
    @test gen.name == "Custom"
    @test gen.nid == 500
    @test length(gen) == 2
    @test total_bp(gen) == 30_000_000
    @test cbp(gen) == UInt32[10_000_000, 30_000_000]
    @test gen.M == 50_000_000

    # one concrete type, ends stored, genetic lengths, optional population size
    @test isconcretetype(Species) && eltype([cattle, sheep, gen]) === Species
    @test cattle.ends === cbp(cattle)
    @test map_length(gen) == [0.2, 0.4]
    @test map_length(cattle) == cattle.chromosome ./ 80_000_000
    micro = GenericSpecies("Micro", [50_000_000, 5_000_000]; morgan = [0.5, 0.5])
    @test map_length(micro) == [0.5, 0.5] && micro.nid == 0
    @test_throws ArgumentError nid(micro)
    @test Cattle().nid == 0 && Cattle().chromosome == cattle.chromosome
    @test map_length(Chicken(; morgan = max.(map_length(Chicken()), 0.5))) |> minimum == 0.5
    @test_throws ArgumentError GenericSpecies("Bad", [1000, 1000]; morgan = [0.1])
    @test_throws ArgumentError GenericSpecies("Bad", [1000, 1000]; morgan = [0.1, 0.0])
    @test_throws ArgumentError GenericSpecies("Big", fill(2_000_000_000, 3))
    @test occursin("not set", sprint(show, MIME("text/plain"), micro))
    old = GenericSpecies("Old", Int32(8), UInt32[100_000, 200_000], UInt32(50_000))
    @test old.nid == 8 && old.M == 50_000 && map_length(old) == [2.0, 4.0]

    @test_throws ArgumentError GenericSpecies("Bad", -1, [1000])
    @test_throws ArgumentError GenericSpecies("Bad", 10, Int[])
    @test_throws ArgumentError GenericSpecies("Bad", 10, [0])
    @test_throws ArgumentError GenericSpecies("Bad", 10, [1000]; M=0)

    for constructor in (Cattle, Sheep, Chicken, Pig, Goat, Horse, Rabbit, Cat, Dog)
        @test_throws ArgumentError constructor(-1)
        @test_throws ArgumentError constructor(100; M=0)
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
    @test_throws ArgumentError Trait("OneCategory", [1.0])
    @test_throws ArgumentError Trait("InfiniteWeight", [Inf, 1.0])
    @test_throws ArgumentError Trait("InfiniteAge"; age=Inf)
    @test_throws ArgumentError Trait("InfiniteMean"; μ=Inf)
end

@testset "Container display" begin
    @test sprint(show, Haplotype(10, 2)) == "Haplotype with 10 loci and 2 haplotypes"
    @test sprint(show, Genotype(10, 2)) == "Genotype with 10 individuals and 2 alleles"
    # Compact forms stay on one line, so they print inside containers.
    for x in (Cattle(10), SNPSet("chip", 5), Trait("milk"), Trait("mast", [0.9, 0.1]))
        @test !occursin('\n', sprint(show, x))
        @test occursin('\n', sprint(show, MIME("text/plain"), x))
    end
    @test 3 in LocusSet("s", [1, 3, 9]) && !(4 in LocusSet("s", [1, 3, 9]))
end

@testset "VariantMap and LocusSet" begin
    chrs = Int8[1, 1, 1, 2, 2]
    pos = UInt32[100, 200, 300, 1150, 1250]   # genome-wide: chromosome 2 starts after 1000
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
    @test sub_vm.pos == UInt32[100, 300, 1250]

    # Coordinate contract: genome-wide, strictly increasing positions
    sp2 = GenericSpecies("Toy", 10, [1000, 1000])
    @test check_coordinates(vm, sp2)
    @test chrpos(vm, sp2) == UInt32[100, 200, 300, 150, 250]
    @test chrpos(vm, sp2, 4) == 150
    within = VariantMap(chrs, [100, 200, 300, 150, 250], sp2; ref = ref, alt = alt, frq = frq)
    @test within.pos == vm.pos && within.chr == vm.chr && within.frq == vm.frq
    derived = VariantMap(pos, sp2)
    @test derived.chr == chrs
    # within-chromosome positions passed as genome-wide ones are rejected
    @test_throws ArgumentError VariantMap(chrs, UInt32[100, 200, 300, 150, 250], ref, alt)
    @test_throws ArgumentError VariantMap(Int8[1, 1], UInt32[100, 100], ['A', 'A'], ['C', 'C'])
    @test_throws ArgumentError VariantMap(Int8[2, 1], UInt32[100, 200], ['A', 'A'], ['C', 'C'])
    @test_throws ArgumentError vm[[3, 1]]
    # a single-chromosome label on a two-chromosome genome fails the species check
    @test_throws ArgumentError check_coordinates(
        VariantMap(fill(Int8(1), 5), pos, ref, alt), sp2)
    @test_throws ArgumentError VariantMap([1, 3], [10, 10], sp2)
    @test_throws ArgumentError VariantMap([1, 2], [10, 2000], sp2)
    @test_throws ArgumentError VariantMap(UInt32[5, 2001], sp2)
    # genetic positions
    gm = VariantMap(chrs, pos, ref, alt; morgan = [0.0, 0.1, 0.2, 0.05, 0.1])
    @test gm[[2, 4]].morgan == [0.1, 0.05]
    @test isempty(vm.morgan)
    @test_throws ArgumentError VariantMap(chrs, pos, ref, alt; morgan = [0.0, 0.3, 0.2, 0.05, 0.1])
    @test_throws ArgumentError VariantMap(chrs, pos, ref, alt; morgan = [0.0, 0.1])
    @test occursin("genome-wide", sprint(show, MIME("text/plain"), gm))
    @test_throws ArgumentError cbp(GenericSpecies("Big", 2, fill(2_000_000_000, 3)))

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
    tbv2 = genotypic_value(hps, qtl2)
    # Ind 1: locus2 (dosage 1) -> 1*2.0 + 0.5 = 2.5; locus3 (dosage 1) -> 1*3.0 = 3.0. Total = 5.5
    # Ind 2: locus2 (dosage 1) -> 1*2.0 + 0.5 = 2.5; locus3 (dosage 1) -> 1*3.0 = 3.0. Total = 5.5
    @test tbv2 ≈ [5.5, 5.5]
    @test genotypic_value(gt, qtl2) ≈ [5.5, 5.5]
    # with dominance, the sum is not a breeding value and tbv refuses it
    @test_throws ArgumentError tbv(hps, qtl2)
    @test genotypic_value(hps, qtl1) == tbv1

    # MultiTraitQTL combining Trait 1 and Trait 2
    mqtl = MultiTraitQTL([qtl1, qtl2])
    @test length(mqtl.trait_names) == 2
    @test mqtl.loci == Int32[1, 2, 3] # union of QTL loci
    @test size(mqtl.additive) == (3, 2)

    @test_throws ArgumentError tbv(hps, mqtl)
    tbv_multi = genotypic_value(hps, mqtl)
    @test size(tbv_multi) == (2, 2)
    @test tbv_multi[:, 1] ≈ tbv1
    @test tbv_multi[:, 2] ≈ tbv2

    @test genotypic_value(gt, mqtl) ≈ tbv_multi
    @test tbv(hps, MultiTraitQTL([qtl1, qtl1])) ≈ [tbv1 tbv1]
    @test_throws ArgumentError TraitQTL("Duplicate", [1, 1], [0.1, 0.2])
end
@testset "Intercept, genotypic value and breeding value" begin
    # One locus in exact Hardy–Weinberg proportions at p = 1/4: 9 : 6 : 1 among 16
    z = [fill(0, 9); fill(1, 6); fill(2, 1)]
    hps = Haplotype(3, 32)
    for (i, zi) in enumerate(z)
        zi >= 1 && (hps[2, 2i-1] = true)
        zi == 2 && (hps[2, 2i] = true)
    end
    a, d, μ, p = 2.0, 1.0, 10.0, 0.25
    qtl = TraitQTL("T", [2], [a]; dominance = [d], intercept = μ)
    G = genotypic_value(hps, qtl)
    @test G == μ .+ z .* a .+ (z .== 1) .* d
    A = breeding_value(hps, qtl)
    α = a + d * (1 - 2p)
    @test A ≈ (z .- 2p) .* α
    @test abs(sum(A)) < 1e-12
    # α is the regression of genotypic value on allele count, and the
    # dominance deviations are uncorrelated with the breeding values
    zc = z .- 2p
    @test sum(zc .* G) / sum(abs2, zc) ≈ α
    D = G .- sum(G) / 16 .- A
    @test abs(sum(A .* D)) < 1e-12
    @test sum(abs2, A) / 16 ≈ 2p * (1 - p) * α^2          # V_A
    @test sum(abs2, D) / 16 ≈ (2p * (1 - p) * d)^2        # V_D
    # frequencies of another reference population
    @test breeding_value(hps, qtl; p = [0.5]) ≈ (z .- 1.0) .* a
    @test_throws DimensionMismatch breeding_value(hps, qtl; p = [0.5, 0.5])
    @test_throws ArgumentError breeding_value(hps, qtl; p = [1.5])
    # additive architecture: breeding value is the centred tbv; the intercept shifts tbv only
    add = TraitQTL("T", [2], [a]; intercept = μ)
    @test tbv(hps, add) == μ .+ z .* a
    @test breeding_value(hps, add) ≈ tbv(hps, add) .- sum(tbv(hps, add)) / 16
    @test breeding_value(hap2id(hps), qtl) ≈ A
    # multi-trait: columns equal the single-trait results, intercepts carried over
    other = TraitQTL("U", [2, 3], [1.0, 5.0]; intercept = -1.0)
    m = MultiTraitQTL([qtl, other])
    @test m.intercept == [μ, -1.0]
    @test genotypic_value(hps, m)[:, 1] ≈ G
    @test genotypic_value(hps, m)[:, 2] ≈ genotypic_value(hps, other)
    @test breeding_value(hps, m)[:, 1] ≈ A
    @test breeding_value(hps, m)[:, 2] ≈ breeding_value(hps, other)
    @test_throws ArgumentError genotypic_value(hps, TraitQTL("far", [4], [1.0]))
    @test_throws ArgumentError TraitQTL("T", [1], [1.0]; intercept = Inf)
    @test_throws ArgumentError MultiTraitQTL(["a", "b"], [1], [1.0 2.0]; intercept = [0.0])
end

@testset "Subsets, dosages and trait columns" begin
    rng = Random.Xoshiro(3)
    bits = rand(rng, Bool, 70, 12)
    hap = Haplotype(BitMatrix(bits))
    sub = select_individuals(hap, [5, 2])
    @test size(sub) == (70, 4)
    @test Matrix(sub) == bits[:, [9, 10, 3, 4]]
    lc = select_loci(hap, LocusSet("s", [3, 66, 70]))
    @test Matrix(lc) == bits[[3, 66, 70], :]
    @test select_loci(hap, 2:4) == bits[2:4, :]
    @test_throws ArgumentError select_individuals(hap, [7])
    @test_throws ArgumentError select_loci(hap, [71])
    @test dosage(hap, 66) == [bits[66, 2i-1] + bits[66, 2i] for i = 1:6]
    Z = dosage(hap, [1, 70])
    @test Z isa Matrix{Int8} && size(Z) == (2, 6)
    @test Z[2, :] == dosage(hap, 70)
    @test dosage(hap2id(hap), [1, 70]) == Z
    @test_throws ArgumentError dosage(hap, 0)
    @test packed(hap) === hap.gt
    @test size(packed(hap), 1) == 128 && !any(packed(hap)[71:end, :])
    cols = trait_columns(Trait("milk"))
    @test cols == (tbv = "tbv_milk", phenotype = "ft_milk", ebv = "ebv_milk")
    @test trait_columns(Trait("mast", [0.9, 0.1])) == trait_columns("mast")
end

@testset "BNGH1 files" begin
    rng = Random.Xoshiro(4)
    bits = rand(rng, Bool, 37, 10)
    hap = Haplotype(BitMatrix(bits))
    pos = UInt32.(cumsum(rand(rng, 1:50, 37)))
    mktempdir() do dir
        f = write_haplotypes(joinpath(dir, "a.bngh"), hap, pos)
        h2, p2 = read_haplotypes(f)
        @test h2 == hap && p2 == pos
        @test packed(h2) == packed(hap)
        vm = VariantMap(pos, GenericSpecies("Toy", 5, [1000, 1000]))
        @test read(write_haplotypes(joinpath(dir, "b.bngh"), hap, vm)) == read(f)
        # the documented layout, byte by byte
        raw = read(f)
        @test raw[1:8] == UInt8['B', 'N', 'G', 'H', '1', 0, 0, 0]
        @test reinterpret(Int64, raw[9:24]) == [37, 10]
        @test reinterpret(UInt32, raw[25:(24+4*37)]) == pos
        @test length(raw) == 24 + 4 * 37 + 37 * 2
        rec = reshape(raw[(25+4*37):end], 2, 37)
        @test all(((rec[(h-1)>>>3+1, l] >>> ((h - 1) & 7)) & 0x01 == 0x01) == bits[l, h]
                  for l = 1:37, h = 1:10)
        @test_throws DimensionMismatch write_haplotypes(joinpath(dir, "c.bngh"), hap, pos[1:3])
        write(joinpath(dir, "bad.bngh"), "not a haplotype file at all, longer than a header")
        @test_throws ArgumentError read_haplotypes(joinpath(dir, "bad.bngh"))
    end
end
