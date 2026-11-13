"""
    struct Species

The genome of a species as the simulation sees it: chromosome lengths and
their genetic lengths. One concrete type serves every species; `Cattle`,
`Sheep`, `Chicken`, `Pig`, `Goat`, `Horse`, `Rabbit`, `Cat` and `Dog` are
constructors with curated autosomal lengths, and [`GenericSpecies`](@ref)
builds a custom genome.

# Fields
- `name::String`: Species or assembly code, e.g. `"BosTau"`.
- `nid::Int32`: Population size, or `0` if not set. **Deprecated**: a
  population size is a property of a population, not of a species. The
  field stays in 0.3 because FisherWright and BnGsim read it, and is to be
  removed once they take the size as their own argument. Prefer the
  constructors without a size.
- `chromosome::Vector{UInt32}`: Chromosome lengths in base pairs (autosomes
  only; sex chromosomes and the mitochondrion are not modelled).
- `M::UInt32`: Base pairs per Morgan, the genome-wide default
  recombination scale (`100_000_000` is 1 cM/Mb).
- `ends::Vector{UInt32}`: Cumulative chromosome ends, the genome-wide
  coordinate of the last base of each chromosome ([`cbp`](@ref)).
- `morgan::Vector{Float64}`: Genetic length of each chromosome in Morgans.
  Defaults to `chromosome ./ M`; give it to model chromosomes whose
  recombination rate per base differs (microchromosomes, or a minimum of
  about 0.5 M per chromosome from the obligate crossover).

The total genome must fit the `UInt32` genome-wide coordinate
(< 4 294 967 296 bp); the constructor checks this.

# Constructor
    Species(name, nid, chromosome; M = 100_000_000, morgan = nothing)

`nid` is a positive integer or `nothing`.

See also [`map_length`](@ref), [`total_bp`](@ref), [`cbp`](@ref).
"""
struct Species
    name::String
    nid::Int32
    chromosome::Vector{UInt32}
    M::UInt32
    ends::Vector{UInt32}
    morgan::Vector{Float64}

    function Species(
        name::AbstractString,
        nid::Union{Nothing,Integer},
        chromosome::AbstractVector{<:Integer};
        M::Integer = 100_000_000,
        morgan::Union{Nothing,AbstractVector{<:Real}} = nothing,
    )
        nid === nothing || nid > 0 || throw(ArgumentError("nid must be positive"))
        nid === nothing || nid <= typemax(Int32) ||
            throw(ArgumentError("nid must not exceed $(typemax(Int32))"))
        isempty(chromosome) && throw(ArgumentError("at least one chromosome length is required"))
        all(>(0), chromosome) || throw(ArgumentError("chromosome lengths must be positive"))
        0 < M <= typemax(UInt32) || throw(ArgumentError("M must be positive and fit UInt32"))
        sum(UInt64, chromosome) <= typemax(UInt32) || throw(ArgumentError(
            "the genome of $name exceeds $(typemax(UInt32)) bp, the range of the UInt32 coordinate"))
        len = convert(Vector{UInt32}, chromosome)
        gl = if morgan === nothing
            Float64.(len) ./ Float64(M)
        else
            length(morgan) == length(len) ||
                throw(ArgumentError("morgan needs one genetic length per chromosome ($(length(len)))"))
            all(m -> isfinite(m) && m > 0, morgan) ||
                throw(ArgumentError("genetic lengths must be finite and positive"))
            convert(Vector{Float64}, morgan)
        end
        new(String(name), Int32(something(nid, 0)), len, UInt32(M), cumsum(len), gl)
    end
end

"""
    GenericSpecies(name, nid, chromosome; M = 100_000_000, morgan = nothing) -> Species
    GenericSpecies(name, chromosome; M = 100_000_000, morgan = nothing) -> Species
    GenericSpecies(name, nid, chromosome, M) -> Species

A [`Species`](@ref) with user-supplied chromosome lengths (bp). The second
form leaves the population size unset; the third is the positional form
that `GenericSpecies` had as a struct before 0.3.

# Examples
```julia
using BnGStructs

sp = GenericSpecies("ModelOrganism", [15_000_000, 25_000_000])
total_bp(sp)     # 40000000
length(sp)       # 2
map_length(sp)   # [0.15, 0.25]

# two chromosomes of very different size, each with 0.5 M
micro = GenericSpecies("Micro", [50_000_000, 5_000_000]; morgan = [0.5, 0.5])
```
"""
GenericSpecies(name::AbstractString, nid::Integer, chromosome::AbstractVector{<:Integer}; kw...) =
    Species(name, nid, chromosome; kw...)
GenericSpecies(name::AbstractString, chromosome::AbstractVector{<:Integer}; kw...) =
    Species(name, nothing, chromosome; kw...)
# The positional form of the struct constructor that `GenericSpecies` had as a
# type before 0.3; FisherWright's examples and tests use it.
GenericSpecies(name::AbstractString, nid::Integer, chromosome::AbstractVector{<:Integer}, M::Integer) =
    Species(name, nid, chromosome; M = M)

"""
    Cattle([nid]; M = 100_000_000, morgan = nothing) -> Species

Create the `Species` of Cattle (*Bos taurus*) genome. `nid`, a population size, is optional and deprecated
(see [`Species`](@ref)).

Autosomal chromosome lengths (29 autosomes) are sourced from NCBI assembly
[GCF_000003055.6](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000003055.6/) (ARS-UCD1.2).

# Arguments
- `nid::Integer`: Population size (optional; must be positive if given).

# Keywords
- `M::Integer`: Base pairs per Morgan (default: `100_000_000`, corresponding to 1 cM/Mb).
- `morgan`: Genetic length of each chromosome in Morgans (default: `chromosome ./ M`).

# Returns
- `Species`: The cattle genome.
"""
function Cattle(nid::Union{Nothing,Integer} = nothing; M::Integer = 100_000_000, morgan = nothing)
    clng = [
        158_534_110,
        136_231_102,
        121_005_158,
        120_000_601,
        120_089_316,
        117_806_340,
        110_682_743,
        113_319_770,
        105_454_467,
        103_308_737,
        106_982_474,
        87_216_183,
        83_472_345,
        82_403_003,
        85_007_780,
        81_013_979,
        73_167_244,
        65_820_629,
        63_449_741,
        71_974_595,
        69_862_954,
        60_773_035,
        52_498_615,
        62_317_253,
        42_350_435,
        51_992_305,
        45_612_108,
        45_940_150,
        51_098_607,
        # 148_823_899, #X
    ]
    return Species("BosTau", nid, clng; M = M, morgan = morgan)
end

"""
    Sheep([nid]; M = 100_000_000, morgan = nothing) -> Species

Create the `Species` of Sheep (*Ovis aries*) genome. `nid`, a population size, is optional and deprecated
(see [`Species`](@ref)).

Autosomal chromosome lengths (26 autosomes) are sourced from NCBI assembly
[GCF_000298735.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000298735.2/) (Oar_v4.0).

# Arguments
- `nid::Integer`: Population size (optional; must be positive if given).

# Keywords
- `M::Integer`: Base pairs per Morgan (default: `100_000_000`, corresponding to 1 cM/Mb).
- `morgan`: Genetic length of each chromosome in Morgans (default: `chromosome ./ M`).

# Returns
- `Species`: The sheep genome.
"""
function Sheep(nid::Union{Nothing,Integer} = nothing; M::Integer = 100_000_000, morgan = nothing)
    clng = [
        275_406_953,
        248_966_461,
        223_996_068,
        119_216_639,
        107_836_144,
        116_888_256,
        100_009_711,
        90_615_088,
        94_583_238,
        86_377_204,
        62_170_480,
        79_028_859,
        82_951_069,
        62_568_341,
        80_783_214,
        71_693_149,
        72_251_135,
        68_494_538,
        60_445_663,
        51_049_468,
        49_987_992,
        50_780_147,
        62_282_865,
        41_976_827,
        45_223_504,
        44_047_080,
        # 135_185_801, #X
    ]
    return Species("OviAri", nid, clng; M = M, morgan = morgan)
end

"""
    Chicken([nid]; M = 100_000_000, morgan = nothing) -> Species

Create the `Species` of Chicken (*Gallus gallus*) genome. `nid`, a population size, is optional and deprecated
(see [`Species`](@ref)).

Autosomal chromosome lengths (39 autosomes) are sourced from NCBI assembly
[GCF_016699485.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_016699485.2/) (bGalGal1.mat.broiler.GRCg7b).

# Arguments
- `nid::Integer`: Population size (optional; must be positive if given).

# Keywords
- `M::Integer`: Base pairs per Morgan (default: `100_000_000`, corresponding to 1 cM/Mb).
- `morgan`: Genetic length of each chromosome in Morgans (default: `chromosome ./ M`).

# Returns
- `Species`: The chicken genome.
"""
function Chicken(nid::Union{Nothing,Integer} = nothing; M::Integer = 100_000_000, morgan = nothing)
    clng = [
        196_449_156,
        149_539_284,
        110_642_502,
        90_861_225,
        59_506_338,
        36_220_557,
        36_382_834,
        29_578_256,
        23_733_309,
        20_453_248,
        19_638_187,
        20_119_077,
        17_905_061,
        15_331_188,
        12_703_657,
        2_706_039,
        11_092_391,
        11_623_896,
        10_455_293,
        14_265_659,
        6_970_754,
        4_686_657,
        6_253_421,
        6_478_339,
        3_067_737,
        5_349_051,
        5_228_753,
        5_437_364,
        726_478,
        755_666,
        2_457_334,
        125_424,
        3_839_931,
        3_469_343,
        554_126,
        358_375,
        157_853,
        667_312,
        177_356,
        # 9_109_940, #W
        # 86_044_486, #Z
        # 16_784, #MT
    ]
    return Species("GalGal", nid, clng; M = M, morgan = morgan)
end

"""
    Pig([nid]; M = 100_000_000, morgan = nothing) -> Species

Create the `Species` of Pig (*Sus scrofa*) genome. `nid`, a population size, is optional and deprecated
(see [`Species`](@ref)).

Autosomal chromosome lengths (18 autosomes) are sourced from NCBI assembly
[GCF_000003025.6](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000003025.6/) (Sscrofa11.1).

# Arguments
- `nid::Integer`: Population size (optional; must be positive if given).

# Keywords
- `M::Integer`: Base pairs per Morgan (default: `100_000_000`, corresponding to 1 cM/Mb).
- `morgan`: Genetic length of each chromosome in Morgans (default: `chromosome ./ M`).

# Returns
- `Species`: The pig genome.
"""
function Pig(nid::Union{Nothing,Integer} = nothing; M::Integer = 100_000_000, morgan = nothing)
    clng = [
        274_330_532,
        151_935_994,
        132_848_913,
        130_910_915,
        104_526_007,
        170_843_587,
        121_844_099,
        138_966_237,
        139_512_083,
        69_359_453,
        79_169_978,
        61_602_749,
        208_334_590,
        141_755_446,
        140_412_725,
        79_944_280,
        63_494_081,
        55_982_971,
        #125_939_595, #X
        #43_547_828,  #Y
    ]
    return Species("SusScr", nid, clng; M = M, morgan = morgan)
end

"""
    Goat([nid]; M = 100_000_000, morgan = nothing) -> Species

Create the `Species` of Goat (*Capra hircus*) genome. `nid`, a population size, is optional and deprecated
(see [`Species`](@ref)).

Autosomal chromosome lengths (29 autosomes) are sourced from NCBI assembly
[GCF_001704415.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001704415.2/) (ARS1).

# Arguments
- `nid::Integer`: Population size (optional; must be positive if given).

# Keywords
- `M::Integer`: Base pairs per Morgan (default: `100_000_000`, corresponding to 1 cM/Mb).
- `morgan`: Genetic length of each chromosome in Morgans (default: `chromosome ./ M`).

# Returns
- `Species`: The goat genome.
"""
function Goat(nid::Union{Nothing,Integer} = nothing; M::Integer = 100_000_000, morgan = nothing)
    clng = [
        157_403_528,
        136_510_947,
        120_038_259,
        120_734_966,
        119_020_588,
        117_642_375,
        108_433_636,
        112_672_867,
        91_568_626,
        101_087_560,
        106_225_002,
        87_277_232,
        83_034_183,
        94_672_733,
        81_904_557,
        79_370_172,
        71_137_785,
        67_275_902,
        62_516_450,
        71_784_255,
        69_425_955,
        60_283_066,
        48_866_549,
        62_310_066,
        42_858_509,
        51_421_553,
        44_709_034,
        44_672_302,
        51_332_696,
    ]
    return Species("CapHir", nid, clng; M = M, morgan = morgan)
end

"""
    Horse([nid]; M = 100_000_000, morgan = nothing) -> Species

Create the `Species` of Horse (*Equus caballus*) genome. `nid`, a population size, is optional and deprecated
(see [`Species`](@ref)).

Autosomal chromosome lengths (31 autosomes) are sourced from NCBI assembly
[GCF_002863925.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_002863925.1/) (EquCab3.0).

# Arguments
- `nid::Integer`: Population size (optional; must be positive if given).

# Keywords
- `M::Integer`: Base pairs per Morgan (default: `100_000_000`, corresponding to 1 cM/Mb).
- `morgan`: Genetic length of each chromosome in Morgans (default: `chromosome ./ M`).

# Returns
- `Species`: The horse genome.
"""
function Horse(nid::Union{Nothing,Integer} = nothing; M::Integer = 100_000_000, morgan = nothing)
    clng = [
        188_260_577,
        121_350_024,
        121_351_753,
        109_462_549,
        96_759_418,
        87_230_776,
        100_787_686,
        97_563_019,
        85_793_548,
        85_155_674,
        61_676_917,
        36_992_759,
        43_784_481,
        94_600_235,
        92_851_403,
        88_962_352,
        80_722_430,
        82_641_348,
        62_681_739,
        65_343_332,
        58_984_458,
        50_928_189,
        55_556_184,
        48_288_683,
        40_282_968,
        43_147_642,
        40_254_690,
        47_348_498,
        34_776_120,
        31_395_959,
        26_001_039,
        #128_206_784,
    ]
    return Species("EquCab", nid, clng; M = M, morgan = morgan)
end

"""
    Rabbit([nid]; M = 100_000_000, morgan = nothing) -> Species

Create the `Species` of Rabbit (*Oryctolagus cuniculus*) genome. `nid`, a population size, is optional and deprecated
(see [`Species`](@ref)).

Autosomal chromosome lengths (21 autosomes) are sourced from NCBI assembly
[GCF_000003625.3](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000003625.3/) (OryCun2.0).

# Arguments
- `nid::Integer`: Population size (optional; must be positive if given).

# Keywords
- `M::Integer`: Base pairs per Morgan (default: `100_000_000`, corresponding to 1 cM/Mb).
- `morgan`: Genetic length of each chromosome in Morgans (default: `chromosome ./ M`).

# Returns
- `Species`: The rabbit genome.
"""
function Rabbit(nid::Union{Nothing,Integer} = nothing; M::Integer = 100_000_000, morgan = nothing)
    clng = [
        194_850_757,
        174_332_312,
        155_691_105,
        91_394_100,
        37_992_211,
        27_502_587,
        173_684_459,
        111_795_807,
        116_251_907,
        47_997_241,
        87_554_214,
        155_355_395,
        143_360_832,
        163_896_628,
        109_054_052,
        84_478_945,
        85_008_467,
        69_800_736,
        57_279_966,
        33_191_332,
        15_578_276,
        #111_700_775,
    ]
    return Species("OryCun", nid, clng; M = M, morgan = morgan)
end

"""
    Cat([nid]; M = 100_000_000, morgan = nothing) -> Species

Create the `Species` of Cat (*Felis catus*) genome. `nid`, a population size, is optional and deprecated
(see [`Species`](@ref)).

Autosomal chromosome lengths (18 autosomes) are sourced from NCBI assembly
[GCF_018350175.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_018350175.1/) (F.catus_Fca126_mat1.0).

# Arguments
- `nid::Integer`: Population size (optional; must be positive if given).

# Keywords
- `M::Integer`: Base pairs per Morgan (default: `100_000_000`, corresponding to 1 cM/Mb).
- `morgan`: Genetic length of each chromosome in Morgans (default: `chromosome ./ M`).

# Returns
- `Species`: The cat genome.
"""
function Cat(nid::Union{Nothing,Integer} = nothing; M::Integer = 100_000_000, morgan = nothing)
    clng = [
        239_360_788,
        169_471_080,
        147_794_258, # A1-3
        200_683_619,
        166_112_658,
        118_522_543,
        114_923_439, # B1-4
        145_515_648,
        129_267_133, # C1-2
        148_844_606,
        90_811_962,
        98_716_366,
        97_937_295, # D1-4
        91_446_128,
        59_620_299,
        59_579_006, # E1-3
        111_732_949,
        52_190_560, # F1-2
        #127_149_045, #X
    ]
    return Species("FelCat", nid, clng; M = M, morgan = morgan)
end

"""
    Dog([nid]; M = 100_000_000, morgan = nothing) -> Species

Create the `Species` of Dog (*Canis lupus familiaris*) genome. `nid`, a population size, is optional and deprecated
(see [`Species`](@ref)).

Autosomal chromosome lengths (38 autosomes) are sourced from NCBI assembly
[GCF_011100685.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_011100685.1/) (UU_Cfam_GSD_1.0).

# Arguments
- `nid::Integer`: Population size (optional; must be positive if given).

# Keywords
- `M::Integer`: Base pairs per Morgan (default: `100_000_000`, corresponding to 1 cM/Mb).
- `morgan`: Genetic length of each chromosome in Morgans (default: `chromosome ./ M`).

# Returns
- `Species`: The dog genome.
"""
function Dog(nid::Union{Nothing,Integer} = nothing; M::Integer = 100_000_000, morgan = nothing)
    clng = [
        123_556_469,
        84_979_418,
        92_479_059,
        89_535_178,
        89_562_946,
        78_113_029,
        81_081_596,
        76_405_709,
        61_171_909,
        70_643_054,
        74_805_798,
        72_970_719,
        64_299_765,
        61_112_200,
        64_676_183,
        60_362_399,
        65_088_165,
        56_472_973,
        55_516_201,
        58_627_490,
        51_742_555,
        61_573_679,
        53_134_997,
        48_566_227,
        51_730_745,
        39_257_614,
        46_662_488,
        41_733_330,
        42_517_134,
        40_643_782,
        39_901_454,
        40_225_481,
        32_139_216,
        42_397_973,
        28_051_305,
        31_223_415,
        30_785_915,
        24_803_098,
        #124_992_030, #X
        # 16_728, #MT
    ]
    return Species("CanFam", nid, clng; M = M, morgan = morgan)
end

function Base.show(io::IO, ::MIME"text/plain", c::Species)
    println(io, lpad("Name: ", 17), c.name)
    println(io, lpad("Population size: ", 17), c.nid > 0 ? c.nid : "not set")
    println(io, lpad("Chromosome bp: ", 17))
    for (i, c) in enumerate(c.chromosome)
        i % 5 == 1 && print(io, ' '^15)
        print(io, lpad(c, 11), ',')
        i % 5 == 0 && println(io)
    end
    length(c.chromosome) % 5 != 0 && println(io)
    println(io, lpad("bp/Morgan: ", 17), Float64(c.M))
    print(io, lpad("Map length: ", 17), round(sum(c.morgan); digits = 3), " M")
end

"""
    name(sp::Species) -> String

Return the species name or code (e.g., `"BosTau"`, `"CanFam"`). Not exported
(the name is too generic); call it as `BnGStructs.name(sp)` or read `sp.name`.

# Arguments
- `sp::Species`: Target species instance.

# Returns
- `String`: Species name identifier.

# Examples
```julia
using BnGStructs

sp = Cattle(100)
BnGStructs.name(sp) # "BosTau"; the same as sp.name
```
"""
name(sp::Species) = sp.name

"""
    nid(sp::Species) -> Int

Return the population size stored in the species, or throw if none is set.
Deprecated together with the `nid` field; see [`Species`](@ref).

# Arguments
- `sp::Species`: Target species instance.

# Returns
- `Int`: Population size.

# Examples
```julia
using BnGStructs

sp = Cattle(100)
nid(sp) # 100
```
"""
function nid(sp::Species)
    sp.nid > 0 || throw(ArgumentError("no population size is set for $(sp.name)"))
    return Int(sp.nid)
end

"""
    chromosome(sp::Species) -> Vector{UInt32}

Return a vector containing the autosomal chromosome lengths in base pairs (bp).

# Arguments
- `sp::Species`: Target species instance.

# Returns
- `Vector{UInt32}`: Autosomal chromosome lengths.

# Examples
```julia
using BnGStructs

sp = Cattle(100)
length(chromosome(sp)) # 29
```
"""
chromosome(sp::Species) = sp.chromosome

"""
    total_bp(sp::Species) -> Integer

Return the total autosomal genome length in base pairs (sum of all chromosome lengths).

# Arguments
- `sp::Species`: Target species instance.

# Returns
- `Integer`: Total length in base pairs.

# Examples
```julia
using BnGStructs

sp = Cattle(100)
total_bp(sp) # sum of all 29 autosome lengths
```
"""
total_bp(sp::Species) = Int(sp.ends[end])

"""
    cbp(sp::Species) -> Vector{UInt32}

Return the cumulative chromosome endpoint positions in base pairs (prefix sums of chromosome lengths).

Useful for translating local chromosome coordinates into genome-wide continuous coordinates.

# Arguments
- `sp::Species`: Target species instance.

# Returns
- `Vector{UInt32}`: Cumulative base-pair endpoints.

# Examples
```julia
using BnGStructs

sp = GenericSpecies("Toy", 10, [1_000, 2_000, 3_000])
cbp(sp) # UInt32[1000, 3000, 6000]
```
"""
cbp(sp::Species) = sp.ends

"""
    M(sp::Species) -> UInt32

Return the recombination scaling factor in base pairs per Morgan. Not exported
(the name is too generic); call it as `BnGStructs.M(sp)` or read `sp.M`.

For example, `100_000_000` base pairs per Morgan corresponds to 1 cM per 1 Mb (1% recombination per Mb).

# Arguments
- `sp::Species`: Target species instance.

# Returns
- `UInt32`: Base pairs per Morgan.

# Examples
```julia
using BnGStructs

sp = Cattle(100)
BnGStructs.M(sp) # 0x05f5e100 (100_000_000); the same as sp.M
```
"""
M(sp::Species) = sp.M
Base.length(sp::Species) = length(sp.chromosome)

"""
    map_length(sp::Species) -> Vector{Float64}

Genetic length of each chromosome in Morgans (the `morgan` field):
`chromosome ./ M` unless lengths were given to the constructor. Its sum is
the genome length in Morgans, the expected number of crossovers per
meiosis under a model without interference.
"""
map_length(sp::Species) = sp.morgan

Base.show(io::IO, c::Species) = print(
    io, "Species \"", c.name, "\" with ", length(c.chromosome), " chromosomes, ",
    Int(c.ends[end]), " bp, ", round(sum(c.morgan); digits = 3), " M",
)
