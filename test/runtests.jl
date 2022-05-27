using Eikonal
using Test

@testset "Eikonal.jl" begin
    @testset "FastSweeping" begin
        m, n = 1000, 1200
        fsm = FastSweeping(m, n)
        fsm.v .= 1

        x0, y0 = 500, 550
        init!(fsm, (x0, y0))
        sweep!(fsm)

        x1, y1 = 50, 50
        @test fsm.t[x1, y1] ≈ sqrt((x1-x0)^2 + (y1-y0)^2)  rtol=1e-2

        r = ray(fsm.t, (x1, y1))
        @test last(r) == (x0, y0)
    end

    @testset "FastMarching" begin
        m, n = 1000, 1200
        fmm = FastMarching(m, n)
        fmm.v .= 1

        x0, y0 = 500, 550
        init!(fmm, (x0, y0))
        march!(fmm)

        x1, y1 = 50, 50
        @test fmm.t[x1, y1] ≈ sqrt((x1-x0)^2 + (y1-y0)^2)  rtol=1e-2

        r = ray(fmm.t, (x1, y1))
        @test last(r) == (x0, y0)
    end
end
