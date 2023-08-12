import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Literate
cd(@__DIR__) do
    @info "running code in `README.jl`"
    include("README.jl")

    Literate.markdown("README.jl", pwd(),
                      flavor=Literate.CommonMarkFlavor())

    # Fix image paths for use in the root directory
    contents = read("README.md", String)
    contents = replace(contents, "![](" => "![](docs/readme/")
    write(joinpath("..", "..", "README.md"), contents)

    Literate.notebook("README.jl", pwd())
end

cd(@__DIR__) do
    run(`pandoc -s --mathjax=https://cdnjs.cloudflare.com/ajax/libs/mathjax/3.2.2/es5/tex-mml-chtml.min.js README.md -o README.html`)
end
