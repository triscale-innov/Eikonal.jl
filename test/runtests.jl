using Eikonal
using Test

@testset "Eikonal.jl" begin
    @testset "FastSweeping" begin
        siz = (1000, 1200)
        fsm = FastSweeping(Float64, siz)
        @test sprint(show, fsm) == "FastSweeping solver on a 1000×1200 grid\n"

        fsm.v .= 1

        x0, y0 = 500, 550
        init!(fsm, (x0, y0))
        sweep!(fsm, verbose=true)

        x1, y1 = 50, 50
        @test fsm.t[x1, y1] ≈ sqrt((x1-x0)^2 + (y1-y0)^2)  rtol=1e-2

        r = ray(fsm.t, (x1, y1))
        @test last(r) == (x0, y0)
    end

    @testset "FastSweeping 3D" begin
        siz = (100, 110, 120)
        fsm = FastSweeping(Float64, siz)
        @test sprint(show, fsm) == "FastSweeping solver on a 100×110×120 grid\n"

        fsm.v .= 1

        source = 50, 55, 60
        init!(fsm, source)
        sweep!(fsm, verbose=true)

        pos = 5, 5, 6
        @test fsm.t[pos...] ≈ sqrt(sum((pos.-source).^2))  rtol=5e-2
    end

    @testset "FastMarching" begin
        m, n = 1000, 1200
        fmm = FastMarching(m, n)
        fmm.v .= 1

        x0, y0 = 500, 550
        init!(fmm, (x0, y0))
        march!(fmm, verbose=true)

        x1, y1 = 50, 50
        @test fmm.t[x1, y1] ≈ sqrt((x1-x0)^2 + (y1-y0)^2)  rtol=1e-2

        r = ray(fmm.t, (x1, y1))
        @test last(r) == (x0, y0)
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

        fsm.t .= min.(fsm.t, 4100)
        r = ray(fsm.t, goal)
        @test last(r) == entrance
    end
end
