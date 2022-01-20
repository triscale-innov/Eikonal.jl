using Eikonal
using Documenter

DocMeta.setdocmeta!(Eikonal, :DocTestSetup, :(using Eikonal); recursive=true)

makedocs(;
    modules=[Eikonal],
    authors="François Févotte <francois.fevotte@triscale-innov.com",
    repo="https://github.com/triscale-innov/Eikonal.jl/blob/{commit}{path}#{line}",
    sitename="Eikonal.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://triscale-innov.github.io/Eikonal.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/triscale-innov/Eikonal.jl",
    devbranch="main",
)
