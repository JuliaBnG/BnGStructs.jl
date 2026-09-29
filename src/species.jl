# Contents
# Species
# Cattle
# Chicken
# Pig
# Sheep
"""
    abstract type Species

Abstract supertype representing species metadata for genetic and genomic simulations.

A `Species` encapsulates species-level genomic parameters including:
- Species name / identifier: [`name`](@ref)
- Reference or population size: [`nid`](@ref)
- Autosomal chromosome lengths in base pairs: [`chromosome`](@ref), [`total_bp`](@ref), [`cbp`](@ref)
- Number of chromosomes: `length(sp)`
- Recombination rate scaling in base pairs per Morgan: [`M`](@ref)

# Subtypes
- [`GenericSpecies`](@ref): Custom species defined by user-supplied chromosome lengths.
- Built-in domestic species with curated autosomal chromosome lengths:
  [`Cattle`](@ref), [`Sheep`](@ref), [`Chicken`](@ref), [`Pig`](@ref),
  [`Goat`](@ref), [`Horse`](@ref), [`Rabbit`](@ref), [`Cat`](@ref), [`Dog`](@ref).
"""
abstract type Species end

"""
    struct GenericSpecies <: Species

User-defined species representation with arbitrary chromosome lengths.

# Fields
- `name::String`: Species identifier.
- `nid::Int32`: Population size (must be positive).
- `chromosome::Vector{UInt32}`: Vector of chromosome lengths in base pairs (bp).
- `M::UInt32`: Base pairs per Morgan (default `100_000_000`, corresponding to 1 cM/Mb).

# Constructors
    GenericSpecies(name::AbstractString, nid::Integer, chromosome::AbstractVector{<:Integer}; M::Integer=100_000_000) -> GenericSpecies

# Examples
```julia
using BnGStructs

sp = GenericSpecies("ModelOrganism", 1000, [15_000_000, 25_000_000]; M=100_000_000)
name(sp)       # "ModelOrganism"
total_bp(sp)   # 40000000
length(sp)     # 2
```
"""
struct GenericSpecies <: Species
    name::String
    nid::Int32
    chromosome::Vector{UInt32} # Changed to UInt32
    M::UInt32 # Changed to UInt32
end

"""
    GenericSpecies(name::AbstractString, nid::Integer, chromosome::AbstractVector{<:Integer}; M::Integer=100_000_000) -> GenericSpecies

Construct a custom `GenericSpecies` with a given `name`, population size `nid`,
chromosome lengths `chromosome` (in bp), and base pairs per Morgan `M`.

# Arguments
- `name::AbstractString`: Name or label for the species.
- `nid::Integer`: Population size (must be positive).
- `chromosome::AbstractVector{<:Integer}`: Vector of positive chromosome lengths in base pairs.

# Keywords
- `M::Integer`: Base pairs per Morgan (default: `100_000_000`).

# Returns
- `GenericSpecies`: Initialized species object.
"""
function GenericSpecies(
    name::AbstractString,
    nid::Integer,
    chromosome::AbstractVector{<:Integer};
    M::Integer=100_000_000,
)
    nid ≤ 0 && error("nid must be positive")
    isempty(chromosome) && error("at least one chromosome length is required")
    all(chromosome .> 0) || error("chromosome lengths must be positive")
    M > 0 || error("M must be positive")
    return GenericSpecies(String(name), Int32(nid), UInt32.(chromosome), UInt32(M))
end

"""
    struct Cattle <: Species

Species definition for domestic cattle (*Bos taurus*).

# Fields
- `name::String`: Species identifier (`"BosTau"`).
- `nid::Int32`: Population size.
- `chromosome::Vector{UInt32}`: Autosomal chromosome lengths (29 autosomes) in base pairs.
- `M::UInt32`: Base pairs per Morgan.
"""
struct Cattle <: Species
    name::String
    nid::Int32
    chromosome::Vector{UInt32} # Changed to UInt32
    M::UInt32 # Changed to UInt32
end

"""
    Cattle(nid::Integer; M::Integer = 100_000_000) -> Cattle

Create a `Cattle` (*Bos taurus*) species object with population size `nid`.

Autosomal chromosome lengths (29 autosomes) are sourced from NCBI assembly
[GCF_000003055.6](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000003055.6/) (ARS-UCD1.2).

# Arguments
- `nid::Integer`: Population size (must be positive).

# Keywords
- `M::Integer`: Base pairs per Morgan (default: `100_000_000`, corresponding to 1 cM/Mb).

# Returns
- `Cattle`: Initialized cattle species instance.
"""
function Cattle(nid::Integer; M::Integer = 100_000_000)
    nid ≤ 0 && error("nid must be positive")
    M > 0 || error("M must be positive")
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
    Cattle("BosTau", Int32(nid), UInt32.(clng), UInt32(M))
end

"""
    struct Sheep <: Species

Species definition for domestic sheep (*Ovis aries*).

# Fields
- `name::String`: Species identifier (`"OviAri"`).
- `nid::Int32`: Population size.
- `chromosome::Vector{UInt32}`: Autosomal chromosome lengths (26 autosomes) in base pairs.
- `M::UInt32`: Base pairs per Morgan.
"""
struct Sheep <: Species
    name::String
    nid::Int32
    chromosome::Vector{UInt32} # Changed to UInt32
    M::UInt32 # Changed to UInt32
end

"""
    Sheep(nid::Integer; M::Integer = 100_000_000) -> Sheep

Create a `Sheep` (*Ovis aries*) species object with population size `nid`.

Autosomal chromosome lengths (26 autosomes) are sourced from NCBI assembly
[GCF_000298735.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000298735.2/) (Oar_v4.0).

# Arguments
- `nid::Integer`: Population size (must be positive).

# Keywords
- `M::Integer`: Base pairs per Morgan (default: `100_000_000`, corresponding to 1 cM/Mb).

# Returns
- `Sheep`: Initialized sheep species instance.
"""
function Sheep(nid::Integer; M::Integer = 100_000_000)
    nid ≤ 0 && error("nid must be positive")
    M > 0 || error("M must be positive")
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
    return Sheep("OviAri", Int32(nid), UInt32.(clng), UInt32(M))
end

"""
    struct Chicken <: Species

Species definition for domestic chicken (*Gallus gallus*).

# Fields
- `name::String`: Species identifier (`"GalGal"`).
- `nid::Int32`: Population size.
- `chromosome::Vector{UInt32}`: Autosomal chromosome lengths (39 autosomes) in base pairs.
- `M::UInt32`: Base pairs per Morgan.
"""
struct Chicken <: Species
    name::String
    nid::Int32
    chromosome::Vector{UInt32} # Changed to UInt32
    M::UInt32 # Changed to UInt32
end

"""
    Chicken(nid::Integer; M::Integer = 100_000_000) -> Chicken

Create a `Chicken` (*Gallus gallus*) species object with population size `nid`.

Autosomal chromosome lengths (39 autosomes) are sourced from NCBI assembly
[GCF_016699485.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_016699485.2/) (bGalGal1.mat.broiler.GRCg7b).

# Arguments
- `nid::Integer`: Population size (must be positive).

# Keywords
- `M::Integer`: Base pairs per Morgan (default: `100_000_000`, corresponding to 1 cM/Mb).

# Returns
- `Chicken`: Initialized chicken species instance.
"""
function Chicken(nid::Integer; M::Integer = 100_000_000)
    nid ≤ 0 && error("nid must be positive")
    M > 0 || error("M must be positive")
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
    return Chicken("GalGal", Int32(nid), UInt32.(clng), UInt32(M))
end

"""
    struct Pig <: Species

Species definition for domestic pig (*Sus scrofa*).

# Fields
- `name::String`: Species identifier (`"SusScr"`).
- `nid::Int32`: Population size.
- `chromosome::Vector{UInt32}`: Autosomal chromosome lengths (18 autosomes) in base pairs.
- `M::UInt32`: Base pairs per Morgan.
"""
struct Pig <: Species
    name::String
    nid::Int32
    chromosome::Vector{UInt32} # Changed to UInt32
    M::UInt32 # Changed to UInt32
end

"""
    Pig(nid::Integer; M::Integer = 100_000_000) -> Pig

Create a `Pig` (*Sus scrofa*) species object with population size `nid`.

Autosomal chromosome lengths (18 autosomes) are sourced from NCBI assembly
[GCF_000003025.6](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000003025.6/) (Sscrofa11.1).

# Arguments
- `nid::Integer`: Population size (must be positive).

# Keywords
- `M::Integer`: Base pairs per Morgan (default: `100_000_000`, corresponding to 1 cM/Mb).

# Returns
- `Pig`: Initialized pig species instance.
"""
function Pig(nid::Integer; M::Integer = 100_000_000)
    nid ≤ 0 && error("nid must be positive")
    M > 0 || error("M must be positive")
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
    return Pig("SusScr", Int32(nid), UInt32.(clng), UInt32(M))
end

"""
    struct Goat <: Species

Species definition for domestic goat (*Capra hircus*).

# Fields
- `name::String`: Species identifier (`"CapHir"`).
- `nid::Int32`: Population size.
- `chromosome::Vector{UInt32}`: Autosomal chromosome lengths (29 autosomes) in base pairs.
- `M::UInt32`: Base pairs per Morgan.
"""
struct Goat <: Species
    name::String
    nid::Int32
    chromosome::Vector{UInt32}
    M::UInt32
end

"""
    Goat(nid::Integer; M::Integer = 100_000_000) -> Goat

Create a `Goat` (*Capra hircus*) species object with population size `nid`.

Autosomal chromosome lengths (29 autosomes) are sourced from NCBI assembly
[GCF_001704415.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001704415.2/) (ARS1).

# Arguments
- `nid::Integer`: Population size (must be positive).

# Keywords
- `M::Integer`: Base pairs per Morgan (default: `100_000_000`, corresponding to 1 cM/Mb).

# Returns
- `Goat`: Initialized goat species instance.
"""
function Goat(nid::Integer; M::Integer = 100_000_000)
    nid ≤ 0 && error("nid must be positive")
    M > 0 || error("M must be positive")
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
    return Goat("CapHir", Int32(nid), UInt32.(clng), UInt32(M))
end

"""
    struct Horse <: Species

Species definition for domestic horse (*Equus caballus*).

# Fields
- `name::String`: Species identifier (`"EquCab"`).
- `nid::Int32`: Population size.
- `chromosome::Vector{UInt32}`: Autosomal chromosome lengths (31 autosomes) in base pairs.
- `M::UInt32`: Base pairs per Morgan.
"""
struct Horse <: Species
    name::String
    nid::Int32
    chromosome::Vector{UInt32}
    M::UInt32
end

"""
    Horse(nid::Integer; M::Integer = 100_000_000) -> Horse

Create a `Horse` (*Equus caballus*) species object with population size `nid`.

Autosomal chromosome lengths (31 autosomes) are sourced from NCBI assembly
[GCF_002863925.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_002863925.1/) (EquCab3.0).

# Arguments
- `nid::Integer`: Population size (must be positive).

# Keywords
- `M::Integer`: Base pairs per Morgan (default: `100_000_000`, corresponding to 1 cM/Mb).

# Returns
- `Horse`: Initialized horse species instance.
"""
function Horse(nid::Integer; M::Integer = 100_000_000)
    nid ≤ 0 && error("nid must be positive")
    M > 0 || error("M must be positive")
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
    return Horse("EquCab", Int32(nid), UInt32.(clng), UInt32(M))
end

"""
    struct Rabbit <: Species

Species definition for European rabbit (*Oryctolagus cuniculus*).

# Fields
- `name::String`: Species identifier (`"OryCun"`).
- `nid::Int32`: Population size.
- `chromosome::Vector{UInt32}`: Autosomal chromosome lengths (21 autosomes) in base pairs.
- `M::UInt32`: Base pairs per Morgan.
"""
struct Rabbit <: Species
    name::String
    nid::Int32
    chromosome::Vector{UInt32}
    M::UInt32
end

"""
    Rabbit(nid::Integer; M::Integer = 100_000_000) -> Rabbit

Create a `Rabbit` (*Oryctolagus cuniculus*) species object with population size `nid`.

Autosomal chromosome lengths (21 autosomes) are sourced from NCBI assembly
[GCF_000003625.3](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000003625.3/) (OryCun2.0).

# Arguments
- `nid::Integer`: Population size (must be positive).

# Keywords
- `M::Integer`: Base pairs per Morgan (default: `100_000_000`, corresponding to 1 cM/Mb).

# Returns
- `Rabbit`: Initialized rabbit species instance.
"""
function Rabbit(nid::Integer; M::Integer = 100_000_000)
    nid ≤ 0 && error("nid must be positive")
    M > 0 || error("M must be positive")
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
    return Rabbit("OryCun", Int32(nid), UInt32.(clng), UInt32(M))
end

"""
    struct Cat <: Species

Species definition for domestic cat (*Felis catus*).

# Fields
- `name::String`: Species identifier (`"FelCat"`).
- `nid::Int32`: Population size.
- `chromosome::Vector{UInt32}`: Autosomal chromosome lengths (18 autosomes: A1-3, B1-4, C1-2, D1-4, E1-3, F1-2) in base pairs.
- `M::UInt32`: Base pairs per Morgan.
"""
struct Cat <: Species
    name::String
    nid::Int32
    chromosome::Vector{UInt32}
    M::UInt32
end

"""
    Cat(nid::Integer; M::Integer = 100_000_000) -> Cat

Create a `Cat` (*Felis catus*) species object with population size `nid`.

Autosomal chromosome lengths (18 autosomes) are sourced from NCBI assembly
[GCF_018350175.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_018350175.1/) (F.catus_Fca126_mat1.0).

# Arguments
- `nid::Integer`: Population size (must be positive).

# Keywords
- `M::Integer`: Base pairs per Morgan (default: `100_000_000`, corresponding to 1 cM/Mb).

# Returns
- `Cat`: Initialized cat species instance.
"""
function Cat(nid::Integer; M::Integer = 100_000_000)
    nid ≤ 0 && error("nid must be positive")
    M > 0 || error("M must be positive")
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
    return Cat("FelCat", Int32(nid), UInt32.(clng), UInt32(M))
end

"""
    struct Dog <: Species

Species definition for domestic dog (*Canis lupus familiaris*).

# Fields
- `name::String`: Species identifier (`"CanFam"`).
- `nid::Int32`: Population size.
- `chromosome::Vector{UInt32}`: Autosomal chromosome lengths (38 autosomes) in base pairs.
- `M::UInt32`: Base pairs per Morgan.
"""
struct Dog <: Species
    name::String
    nid::Int32
    chromosome::Vector{UInt32}
    M::UInt32
end

"""
    Dog(nid::Integer; M::Integer = 100_000_000) -> Dog

Create a `Dog` (*Canis lupus familiaris*) species object with population size `nid`.

Autosomal chromosome lengths (38 autosomes) are sourced from NCBI assembly
[GCF_011100685.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_011100685.1/) (UU_Cfam_GSD_1.0).

# Arguments
- `nid::Integer`: Population size (must be positive).

# Keywords
- `M::Integer`: Base pairs per Morgan (default: `100_000_000`, corresponding to 1 cM/Mb).

# Returns
- `Dog`: Initialized dog species instance.
"""
function Dog(nid::Integer; M::Integer = 100_000_000)
    nid ≤ 0 && error("nid must be positive")
    M > 0 || error("M must be positive")
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
    return Dog("CanFam", Int32(nid), UInt32.(clng), UInt32(M))
end

function Base.show(io::IO, c::Species)
    println(io, lpad("Name: ", 17), c.name)
    println(io, lpad("Population size: ", 17), c.nid)
    println(io, lpad("Chromosome bp: ", 17))
    for (i, c) in enumerate(c.chromosome)
        i % 5 == 1 && print(io, ' '^15)
        print(io, lpad(c, 11), ',')
        i % 5 == 0 && println(io)
    end
    length(c.chromosome) % 5 != 0 && println(io)
    print(io, lpad("bp/Morgan: ", 17), Float64(c.M))
end

"""
    name(sp::Species) -> String

Return the species name or code (e.g., `"BosTau"`, `"CanFam"`).

# Arguments
- `sp::Species`: Target species instance.

# Returns
- `String`: Species name identifier.

# Examples
```julia
using BnGStructs

sp = Cattle(100)
name(sp) # "BosTau"
```
"""
name(sp::Species) = sp.name

"""
    nid(sp::Species) -> Int

Return the species population size.

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
nid(sp::Species) = Int(sp.nid)

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
total_bp(sp::Species) = sum(sp.chromosome)

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
cbp(sp::Species) = cumsum(sp.chromosome)

"""
    M(sp::Species) -> UInt32

Return the recombination scaling factor in base pairs per Morgan.

For example, `100_000_000` base pairs per Morgan corresponds to 1 cM per 1 Mb (1% recombination per Mb).

# Arguments
- `sp::Species`: Target species instance.

# Returns
- `UInt32`: Base pairs per Morgan.

# Examples
```julia
using BnGStructs

sp = Cattle(100)
M(sp) # 0x05f5e100 (100_000_000)
```
"""
M(sp::Species) = sp.M
Base.length(sp::Species) = length(sp.chromosome)
