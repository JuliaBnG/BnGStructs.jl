using Documenter
using BnGStructs

# JuliaBnG Sepia and Black themes, added to Documenter's theme picker.
# The files in src/assets/juliabng/ are built by docs-themes/build.jl in the
# juliabng.github.io repository; rebuild there rather than editing them here.
const JULIABNG_ASSETS = [
    [asset("assets/juliabng/$t.css"; class = :css, islocal = true,
           attributes = Dict(Symbol("data-theme-name") => t))
     for t in ("juliabng-sepia", "juliabng-black")];
    "assets/juliabng/themes.js"
]

makedocs(
    modules = [BnGStructs],
    sitename = "BnGStructs.jl",
    authors = "Xijiang Yu",
    format = Documenter.HTML(assets = JULIABNG_ASSETS),
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
