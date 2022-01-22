include("utils.jl")

fm = setup!(FastMarching, TestCase)
fm′ = deepcopy(fm)

@time march!(fm)

plt = background(TestCase, fm, "Fast Marching Method", tmax=fm.t[3050, 475])
rays(TestCase, fm, plt, "fmm")

anim(TestCase, fm′,
     tmax = fm.t[3050, 475])
