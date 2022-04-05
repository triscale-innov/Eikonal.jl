using Eikonal
using Plots
using Printf
using PrettyTables

function _savefig(plt, fname::String)
    mktempdir() do path
        tmp = joinpath(path, fname)
        savefig(plt, tmp)
        mv(tmp, fname, force=true)
    end
end


begin
    struct TestCase end

    function setup!(solver, ::Type{TestCase})
        sol = solver(4000, 3000)
        @info "Initializing solver"  method=solver case="TestCase" npixels=prod(size(sol.t))

        init!(sol, (1000, 1000))
        init!(sol, (3000, 2000))

        sol.v .= 1

        # Slow zone
        sol.v[1300:2000, 1500:2000] .= 5

        # Labyrinth (with "unreachable" walls)
        CI = CartesianIndex
        sol.v[CI(3000,500):CI(3200,505)] .= 10000
        sol.v[CI(3000,400):CI(3005,500)] .= 10000
        sol.v[CI(3000,400):CI(3300,405)] .= 10000
        sol.v[CI(3300,400):CI(3305,600)] .= 10000
        sol.v[CI(2900,600):CI(3305,605)] .= 10000
        sol.v[CI(2900,300):CI(2905,600)] .= 10000

        sol
    end

    function background(::Type{TestCase}, sol, title; tmax=maximum(sol.t))
        plt = plot(title=title, dpi=300)
        contour!(min.(sol.t',tmax), levels=20, fill=true, c=:coolwarm)
    end

    function rays(::Type{TestCase}, sol, plt, output)
        for pos in ((1500, 2990),
                    (1930, 1600),
                    (1900, 1600),
                    (1500, 1950),
                    (3500,  300),
                    (3500,  300),
                    (3050,  475))
            r = ray(sol.t, pos, ρ=1e-2)
            plot!(plt, first.(r), last.(r),
                  linewidth=2, linecolor=:green3, label=nothing)
        end
        _savefig(plt, "$output.svg")
        _savefig(plt, "$output.png")
    end
end



begin
    struct Sandwich end

    function setup!(solver, ::Type{Sandwich})
        # Geometric domain
        m = 2000  # × 10⁻⁴ m  = 20.0 cm
        n = 1150  # × 10⁻⁴ m  = 11.5 cm
        sol = solver(m, n)

        @info "Initializing solver"  method=solver case="Sandwich" npixels=prod(size(sol.t))

        # Boundary condition of the eikonal equation
        # (= source of the wave)
        init!(sol, (1001, 1150));

        # Material speeds [m/s]
        c = (
            steel     = 5900,
            cardboard = 2300,
            copper    = 3600,
            iron      = 5900,
            oil       = 1400,
        )

        # Multiplicative factor on speeds to get arrival times in µs
        units = prod((
            1e-4, # m/cell
            1e6   # µs/s
        ))

        # Inverse speed field
        sol.v .= NaN
        sol.v[:,    1: 200] .= units / c.iron
        sol.v[:,  201: 500] .= units / c.oil
        sol.v[:,  501: 550] .= units / c.cardboard
        sol.v[:,  551: 750] .= units / c.copper
        sol.v[:,  751: 800] .= units / c.cardboard
        sol.v[:,  801:1100] .= units / c.oil
        sol.v[:, 1101:1150] .= units / c.steel

        sol
    end

    function background(::Type{Sandwich}, sol, title; tmax=maximum(sol.t))
        plt = plot(title=title, dpi=300)

        # Level curves of arrival times
        (m, n) = size(sol.v)
        X = collect((0:m) .* 0.1)
        Y = collect((0:n) .* 0.1)
        contour!(X, Y, min.(sol.t', tmax), levels=30, fill=true, c=:coolwarm)

        # Homogenous region boundaries
        for y in (20, 50, 55, 75, 80, 110)
            plot!(plt, [0, 200], [y, y], label=nothing, color=:black, line=:dash)
        end

        plt
    end

    function rays(::Type{Sandwich}, sol, plt, output)
        @info "Retrieving arrival times and shortest paths for sensors"

        # Retrieve arrival times for each sensor
        table = []
        for sensor in (("acier_7",   2000, 1126),
                       # ("huile_6",   2000, 951),
                       ("carton_5",  2000, 771),
                       ("cuivre_4",  2000, 651),
                       ("carton_3",  2000, 526),
                       # ("huile_2",   2000, 351),
                       ("fer_1",     2000, 101),
                       ("bottom",    1001,   2),
                       )
            (name, x, y) = sensor
            push!(table, hcat(name, 0.1*(x-1), 0.1*(y-1), sol.t[x, y]))

            # and reconstruct the rays using a gradient descent algorithm
            r = ray(sol.t, (x,y))
            plot!(plt, 0.1*(first.(r).-1), 0.1*(last.(r).-1),
                  linewidth=2, linecolor=:green3, label=nothing)
        end
        pretty_table(reduce(vcat, table),
                     header = (["Sensor", "x",    "y",    "Arrival time"],
                               ["name",   "[mm]", "[mm]", "[µs]"]),
                     alignment = [:l,:r,:r,:r],
                     formatters = ft_printf("%5.2f", 2:4))

        _savefig(plt, "$output.svg")
        _savefig(plt, "$output.png")
    end
end










function anim(typ, fs::FastSweeping;
              title = "Fast Sweeping Method",
              output = "fsm.gif",
              tmax)
    @info "Generating animation"
    anim = Animation()
    while true
        converged = sweep!(fs, verbose=true, nsweeps=1, epsilon=1e-3)
        (k, l) = fs.iter[]
        l = collect("⬈⬉⬋⬊")[l]
        plt = background(typ, fs, title; tmax=tmax)
        annotate!(plt, 0,0, (" iter $k$l", :bottom, :left))
        frame(anim, plt)
        converged && break
    end
    gif(anim, output, fps=1)
end

function anim(typ, fm::FastMarching;
              title = "FastMarching Method",
              output = "fmm.gif",
              tmax)
    @info "Generating animation"
    a = Animation()
    δt = tmax / 30
    for horizon in δt:δt:tmax
        converged = march!(fm; tmax=horizon, verbose=true)
        plt = background(typ, fm, title; tmax=tmax)
        annotate!(plt, 0,0, (@sprintf(" t = %.2f", horizon), :bottom, :left))
        frame(a, plt)

        converged && break
    end
    gif(a, "fmm.gif", fps=5)
end




run!(fsm::FastSweeping) = sweep!(fsm, verbose=true, epsilon=1e-5)
run!(fmm::FastMarching) = march!(fmm, verbose=true)

# Indice (i,j) du pixel le plus proche de la position (x,y) dans la direction δ
#
# (i,j) = pos2ind((x,y), δ)
# - si δ > 0, le centre du pixel (i,j) est au nord-est  de la position (x,y)
# - si δ < 0, le centre du pixel (i,j) est au sud-ouest de la position (x,y)
pos2ind(sw, pos, δ) = ceil.(Int,  δ .+ (pos .- sw) ./ h)
pos2ind(sw, pos)    = round.(Int,      (pos .- sw) ./ h)

function rasterize(input, h, materials::Dict{Symbol, T}) where {T}
    sw = Tuple(input["bounding_box"]["sw"])
    ne = Tuple(input["bounding_box"]["ne"])

    (Nx, Ny) = pos2ind(sw, ne, -0.5)
    grid = fill(materials[:UNKNOWN], Nx, Ny)

    k = 0
    for rect in input["rects"]
        (i1, j1) = pos2ind(sw, Tuple(rect["sw"]),  0.5)
        (i2, j2) = pos2ind(sw, Tuple(rect["ne"]), -0.5)

        grid[i1:i2, j1:j2] .= materials[Symbol(rect["material"])]
    end

    grid
end
