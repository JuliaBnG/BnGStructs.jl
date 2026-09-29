using Documenter
using BnGStructs

makedocs(
    modules = [BnGStructs],
    sitename = "BnGStructs.jl",
    authors = "Xijiang Yu",
    format = Documenter.HTML(),
    checkdocs = :exports,
    pages = [
        "Home" => "index.md",
        "API Reference" => "api.md",
    ],
)

deploydocs(
    repo = "github.com/JuliaBnG/BnGStructs.jl.git",
    deploy_repo = "github.com/JuliaBnG/juliabng.github.io.git",
    dirname = "BnGStructs",
)
