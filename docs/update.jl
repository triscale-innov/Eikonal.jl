#!/usr/bin/env -S julia --startup-file=no

import Pkg
Pkg.activate(@__DIR__)
Pkg.resolve()

using Literate

function update(fname, target=nothing)
    @info "running code in `$fname`"
    include(abspath(fname))

    Literate.markdown(fname, pwd(),
                      flavor=Literate.CommonMarkFlavor())

    Literate.notebook(fname, pwd())

    name, _ = splitext(fname)
    run(`pandoc -s --mathjax=https://cdnjs.cloudflare.com/ajax/libs/mathjax/3.2.2/es5/tex-mml-chtml.min.js $name.md -o $name.html`)

    if target !== nothing
        @info "writing result to target directory `$target`"
        contents = read("$name.md", String)
        cd(target) do
            write("$name.md", contents)
            run(`pandoc -s --mathjax=https://cdnjs.cloudflare.com/ajax/libs/mathjax/3.2.2/es5/tex-mml-chtml.min.js $name.md -o $name.html`)
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    fname = get(ARGS, 1, "README.jl")
    target = get(ARGS, 2, nothing)
    update(fname, target)
end

true
