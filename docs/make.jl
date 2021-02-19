using OscillatoryIntegralsODE
using Documenter

DocMeta.setdocmeta!(OscillatoryIntegralsODE, :DocTestSetup, :(using OscillatoryIntegralsODE); recursive=true)

makedocs(;
    modules=[OscillatoryIntegralsODE],
    authors="Zack Li",
    repo="https://github.com/xzackli/OscillatoryIntegralsODE.jl/blob/{commit}{path}#{line}",
    sitename="OscillatoryIntegralsODE.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://xzackli.github.io/OscillatoryIntegralsODE.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "API" => "api_index.md"
    ],
)

deploydocs(;
    repo="github.com/xzackli/OscillatoryIntegralsODE.jl",
    devbranch = "main"
)
