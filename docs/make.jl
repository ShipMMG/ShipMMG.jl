using ShipMMG
using Documenter

DocMeta.setdocmeta!(ShipMMG, :DocTestSetup, :(using ShipMMG); recursive = true)

makedocs(;
    modules = [ShipMMG],
    authors = "Taiga MITSUYUKI",
    repo = "https://github.com/ShipMMG/ShipMMG.jl/blob/{commit}{path}#{line}",
    sitename = "ShipMMG.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://ShipMMG.github.io/ShipMMG.jl",
        assets = String[],
    ),
    pages = ["Home" => "index.md"],
)

deploydocs(; repo = "github.com/ShipMMG/ShipMMG.jl")
