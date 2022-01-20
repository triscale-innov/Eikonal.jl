using Eikonal
using Plots
using Printf

m = 4000
n = 3000
fm = FastMarching(m, n)

init!(fm, (1000, 1000))
init!(fm, (3000, 2000))

fm.v .= 1

# Slow zone
fm.v[1300:2000, 1500:2000] .= 5

# Labyrinth (with "unreachable" walls)
const CI = CartesianIndex
fm.v[CI(3000,500):CI(3200,505)] .= 10000
fm.v[CI(3000,400):CI(3005,500)] .= 10000
fm.v[CI(3000,400):CI(3300,405)] .= 10000
fm.v[CI(3300,400):CI(3305,600)] .= 10000
fm.v[CI(2900,600):CI(3305,605)] .= 10000
fm.v[CI(2900,300):CI(2905,600)] .= 10000
fm′ = deepcopy(fm)


# Post-processing & visualization
function plot_contour(fm; tmax=maximum(fm.t))
    plt = plot(title="Test FMM", dpi=300)
    contour!(min.(fm.t', tmax), levels=20, fill=true, c=:coolwarm)
end

@info "Computing arrival times using the FMM"
@time march!(fm)

@info "Generating image"
begin
    plt = plot_contour(fm, tmax=fm.t[3050, 475])
    for pos in ((1500, 2990),
                (1930, 1600),
                (1900, 1600),
                (1500, 1950),
                (3500,  300),
                (3500,  300),
                (3050,  475))
        r = ray(fm.t, pos, ρ=0.01)
        plot!(plt, first.(r), last.(r),
              linewidth=2, linecolor=:green3, label=nothing)
    end

    savefig(plt, "tmp.svg"); mv("tmp.svg", "fmm.svg", force=true)
    savefig(plt, "tmp.png"); mv("tmp.png", "fmm.png", force=true)
end

@info "Generating animation"
begin
    anim = Animation()
    for horizon in 50:50:2600
        converged = march!(fm′; tmax=horizon, verbose=true)
        frame(anim, plot_contour(fm′, tmax=fm.t[3050, 475]))
        converged && break
    end
    gif(anim, "fmm.gif", fps=5)
end
