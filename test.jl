using FastSweeping
using Plots

m = 4000
n = 3000
t = fill(Inf, m, n)
v = similar(t)

t[1000, 1000] = 0
t[3000, 2000] = 0

v .= 1
v[1300:2000, 1500:2000] .= 5

@time sweep!(t, v)

plt = plot(title="Fast Sweeping Method", dpi=300)
contour!(t')

plot_ray!(plt, r) = scatter!(plt, first.(r), last.(r), markersize=0.01, label=nothing)
r = ray(t, (1500, 2990)); plot_ray!(plt, r)
r = ray(t, (1930, 1600)); plot_ray!(plt, r)
r = ray(t, (1900, 1600)); plot_ray!(plt, r)
r = ray(t, (1500, 1950)); plot_ray!(plt, r)
r = ray(t, (3500,  300)); plot_ray!(plt, r)
savefig(plt, "fast_sweeping.svg")
savefig(plt, "fast_sweeping.png")

plt
