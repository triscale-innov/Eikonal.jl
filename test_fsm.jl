using Eikonal
using Plots

m = 4000
n = 3000
t = fill(Inf, m+1, n+1)
v = similar(t, m, n)

t[1000, 1000] = 0
t[3000, 2000] = 0

v .= 1

# Slow zone
v[1300:2000, 1500:2000] .= 5

# Labyrinth (with "unreachable" walls)
const CI = CartesianIndex
v[CI(3000,500):CI(3200,505)] .= 10000
v[CI(3000,400):CI(3005,500)] .= 10000
v[CI(3000,400):CI(3300,405)] .= 10000
v[CI(3300,400):CI(3305,600)] .= 10000
v[CI(2900,600):CI(3305,605)] .= 10000
v[CI(2900,300):CI(2905,600)] .= 10000


@time sweep!(t, v, verbose=true)

plt = plot(title="Fast Sweeping Method", dpi=300)

# rescale t to avoid arbitrarily large values coming from the "labyrinth"
tmax = t[3050, 475]
contour!(min.(t',tmax), levels=20, fill=true, c=:coolwarm)

for pos in ((1500, 2990),
            (1930, 1600),
            (1900, 1600),
            (1500, 1950),
            (3500,  300),
            (3500,  300),
            (3050,  475))
    r = ray(t, pos, ρ=1e-2)
    plot!(plt, first.(r), last.(r),
          linewidth=2, linecolor=:green3, label=nothing)
end

savefig(plt, "tmp.svg"); mv("tmp.svg", "fsm.svg", force=true)
savefig(plt, "tmp.png"); mv("tmp.png", "fsm.png", force=true)
