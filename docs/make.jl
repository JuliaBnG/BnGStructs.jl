pushfirst!(LOAD_PATH, joinpath(@__DIR__, ".."))

using Documenter
using BnGStructs

makedocs(
    modules=[BnGStructs],
    sitename="BnGStructs.jl",
    authors="Xijiang Yu",
    format=Documenter.HTML(prettyurls=false),
    checkdocs=:exports,
    pages=[
        "Home" => "index.md",
        "API Reference" => "api.md",
    ],
)
