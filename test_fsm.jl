include("utils.jl")

fs  = setup!(FastSweeping, TestCase)
fs′ = deepcopy(fs)

@time sweep!(fs, verbose=true)

plt = background(TestCase, fs, "Fast Sweeping Method", tmax=fs.t[3050, 475])
rays(TestCase, fs, plt, "fsm")

anim(TestCase, fs′,
     tmax = fs.t[3050, 475])
