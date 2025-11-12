# BnGStructs.jl

BnGStructs.jl is a Julia package providing efficient data structures for representing and manipulating SNP haplotypes and genotypes in population genetics and breeding simulations. It is designed for use in genomic selection, simulation, and quantitative genetics workflows.

## Features
- `Haplotype` struct: compact storage of haplotype-major SNP data using `BitMatrix`
- `Genotype` struct: compact storage of ID-major SNP data using `BitMatrix`
- Conversion utilities between haplotype-major and ID-major representations
- Type-stable, memory-efficient, and compatible with multi-threaded operations

## Installation
Until registered, install via path or Git URL:
```julia
] #to enter pkg environment
add BnGStructs
```

## Usage Example
```julia
using BnGStructs

# Create a Haplotype with 1000 loci and 200 haplotypes (must be even)
hap = Haplotype(1000, 200)

# Convert to Genotype (requires mtGEBV dependency)
gt = hap2id(hap)

# Access underlying BitMatrix
gtmat = gt.gt

# Define a Cattle population
cattle = Cattle(1000)
```

## Requirements for Registration
- [x] MIT License (`LICENSE` file present)
- [x] Julia-specific `.gitignore`
- [x] `Project.toml` with name, uuid, version, authors, [compat] with `julia = "1"` or similar
- [x] All source in `src/`, entrypoint `BnGStructs.jl`
- [x] Minimal tests in `test/runtests.jl` (recommended for registration)
- [x] No `Manifest.toml` committed (recommended for libraries)
- [x] Public repository (e.g., https://github.com/JuliaBnG/BnGStructs.jl)

## How to Register
1. Tag a release (e.g., `v0.1.0`) and push to GitHub
2. Install the [Registrator GitHub App](https://github.com/JuliaRegistries/Registrator.jl#github-app) on your repo
3. Comment `@JuliaRegistrator register` on the tagged commit or open an issue
4. Address any registry review feedback (compat, tests, etc.)

## License
MIT. See `LICENSE`.
