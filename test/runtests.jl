using Eikonal
using Test

using Eikonal: NearestMin

dist(x, y) = maximum(x.-y)

@testset "Eikonal.jl" begin
    @testset "FastSweeping" begin
        siz = (1100, 1200)
        fsm = FastSweeping(Float64, siz)
        @test sprint(show, fsm) == "FastSweeping solver on a 1100×1200 grid"

        fsm.v .= 1

        source = 500, 550
        init!(fsm, source)
        sweep!(fsm, verbose=true)

        pos = 50, 50
        @test fsm.t[pos...] ≈ sqrt(sum((pos.-source).^2))  rtol=1e-2

        r = ray(fsm.t, pos)
        @test dist(last(r), source) <= 1

        r = ray(fsm.t, pos, NearestMin())
        @test dist(last(r), source) <= 1
    end

    @testset "FastSweeping 3D" begin
        siz = (100, 110, 120)
        fsm = FastSweeping(Float64, siz)
        @test sprint(show, fsm) == "FastSweeping solver on a 100×110×120 grid"

        fsm.v .= 1

        source = 50, 55, 60
        init!(fsm, source)
        sweep!(fsm, verbose=true)

        pos = 5, 5, 6
        @test fsm.t[pos...] ≈ sqrt(sum((pos.-source).^2))  rtol=5e-2

        r = ray(fsm.t, pos)
        @test dist(last(r), source) <= 1

        r = ray(fsm.t, pos, NearestMin())
        @test dist(last(r), source) <= 1
    end

    @testset "FastMarching" begin
        m, n = 1000, 1200
        fmm = FastMarching(m, n)
        fmm.v .= 1

        source = 500, 550
        init!(fmm, source)
        march!(fmm, verbose=true)

        pos = 50, 50
        @test fmm.t[pos...] ≈ sqrt(sum((pos.-source).^2))  rtol=1e-2

        r = ray(fmm.t, pos)
        @test dist(last(r), source) <= 1
    end

    @testset "Utilities" begin
        img = joinpath(@__DIR__, "..", "docs", "readme", "maze.png")
        fsm = FastSweeping(img, ["black"=>Inf, "white"=>1.0])

        entrance = (10, 100)
        init!(fsm, entrance)
        sweep!(fsm, epsilon=1e-5)

        goal = size(fsm.v) .- (10, 100)
        t = vertex2cell(fsm.t)
        @test fsm.t[goal...] ≈ 3937        rtol=1e-2
        @test fsm.t[goal...] ≈ t[goal...]  rtol=1e-2

        r = ray(fsm.t, goal)
        @test dist(last(r), entrance) <= 1

        r = ray(fsm.t, goal, NearestMin())
        @test dist(last(r), entrance) <= 1
    end
end
