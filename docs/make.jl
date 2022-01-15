using FastSweeping
using Documenter

DocMeta.setdocmeta!(FastSweeping, :DocTestSetup, :(using FastSweeping); recursive=true)

makedocs(;
    modules=[FastSweeping],
    authors="François Févotte <francois.fevotte@triscale-innov.com",
    repo="https://github.com/triscale-innov/FastSweeping.jl/blob/{commit}{path}#{line}",
    sitename="FastSweeping.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://triscale-innov.github.io/FastSweeping.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/triscale-innov/FastSweeping.jl",
    devbranch="main",
)
