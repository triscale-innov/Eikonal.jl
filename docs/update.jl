import Pkg
Pkg.activate(@__DIR__)
Pkg.resolve()

using Literate

function update(fname)
    @info "running code in `$fname`"
    include(fname)

    Literate.markdown(fname, pwd(),
                      flavor=Literate.CommonMarkFlavor())

    Literate.notebook(fname, pwd())

    name, _ = splitext(fname)
    run(`pandoc -s --mathjax=https://cdnjs.cloudflare.com/ajax/libs/mathjax/3.2.2/es5/tex-mml-chtml.min.js $name.md -o $name.html`)
end

true
