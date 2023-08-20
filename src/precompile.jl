using PrecompileTools

@setup_workload let
    σ₂ = rand(10, 10)
    σ₃ = rand(10, 10, 10)
    @compile_workload begin
        fs₂ = FastSweeping(Float64, (10, 10))
        fs₂ = FastSweeping(σ₂)
        init!(fs₂, (5, 5))
        sweep!(fs₂, verbose=false, epsilon=1e-5, nsweeps=1)
        sweep!(fs₂)
        ray(fs₂.t, (6, 6))

        fs₃ = FastSweeping(Float64, (10, 10, 10))
        fs₃ = FastSweeping(σ₃)
        init!(fs₃, (5, 5, 5))
        sweep!(fs₃, verbose=false, epsilon=1e-5, nsweeps=1)
        sweep!(fs₃)
        ray(fs₃.t, (6, 6, 6))

        fm₂ = FastMarching(10, 10)
        fm₂.v .= 1
        init!(fm₂, (5, 5))
        march!(fm₂, verbose=false, tmax=1e6, itermax=10)
        march!(fm₂)
    end
end
