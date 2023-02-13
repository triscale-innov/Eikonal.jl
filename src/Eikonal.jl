module Eikonal

using DataStructures
using Images
using LinearAlgebra
using Printf

export FastSweeping, sweep!
export FastMarching, init!, march!
export ray
export vertex2cell


struct FastSweeping{T}
    t      :: Matrix{T}
    v      :: Matrix{T}
    change :: Vector{T}
    iter   :: Base.RefValue{Tuple{Int, Int}}  # (k, l)
end

function Base.show(io::IO, fs::FastSweeping)
    (m,n) = size(fs.v)
    println(io, "FastSweeping solver on a $m×$n grid")
end


function FastSweeping(m::Int, n::Int, T=Float64)
    t = fill(typemax(T), m+1, n+1)
    v = fill(typemax(T), m,   n)
    change = ones(T, 4)
    iter = Ref((1, 0))
    FastSweeping(t, v, change, iter)
end

function init!(fs :: FastSweeping, source)
    (i, j) = source
    fs.t[i, j] = 0
end

# Fast Sweeping Method
function sweep!(fs :: FastSweeping{T};
                verbose = false,
                epsilon = eps(T),
                nsweeps = typemax(Int),
                ) where {T}
    (t, v, change, iter) = (fs.t, fs.v, fs.change, fs.iter)
    (m, n) = size(v)

    ordered_indices(N, δ) = δ>0 ? (2:1:N+1) : (N:-1:1)
    quadrants = ((1,1), (-1,1), (-1,-1), (1, -1))

    k = first(iter[]) # Global iter number
    l = last(iter[])  # Quadrant index
    sweep = 0         # Sweep number
    @inbounds while true
        l += 1
        sweep += 1
        if l == 5
            l = 1
            k += 1
        end
        iter[] = (k, l)

        (δᵢ,δⱼ) = quadrants[l]
        maxchange = change[l] = zero(T)
        maximum(change) <= epsilon && return true

        sᵢ = δᵢ<0 ? 0 : 1   # Offset between vertex-based indices (in `t`)
        sⱼ = δⱼ<0 ? 0 : 1   # and cell-based indices (in `v`)

        for j in ordered_indices(n, δⱼ)
            for i in ordered_indices(m, δᵢ)
                vᵢⱼ = v[i-sᵢ, j-sⱼ]
                tᵢ  = t[i-δᵢ, j]
                tⱼ  = t[i, j-δⱼ]

                # Hyp: sign(∇x) == sign(δi) && sign(∇y) == sign(δj)
                #
                # a⋅tᵢⱼ² + b⋅tᵢⱼ + c = 0
                a = 2
                b = -2*(tᵢ + tⱼ)
                c = tᵢ^2 + tⱼ^2 - vᵢⱼ^2
                Δ = b^2 - 4*a*c

                t₁ = t₂ = typemax(T)
                if Δ>0                         # TODO: really need to compute both roots?
                    t₁ = (-b + √(Δ)) / 2a
                    t₁ > max(tᵢ, tⱼ) || (t₁ = typemax(T))

                    t₂ = (-b - √(Δ)) / 2a
                    t₂ > max(tᵢ, tⱼ) || (t₂ = typemax(T))
                end

                # Hyp: ∇y == 0  or  ∇x == 0
                t₃ = vᵢⱼ + tᵢ
                t₄ = vᵢⱼ + tⱼ

                # Update if necessary
                tmin = min(t₁, t₂, t₃, t₄)
                if tmin < t[i,j]
                    maxchange = max(maxchange, (t[i,j]-tmin)/tmin)
                    t[i,j] = tmin
                end
            end

        end
        verbose && println("iter $k, sweep $l: change = $maxchange")
        change[l] = maxchange

        sweep >= nsweeps && return false
    end

    error("Max number of iterations reached!")
end





function update(t::AbstractMatrix{T}, v, (i,j), (δᵢ, δⱼ)) where {T}
    sᵢ = δᵢ<0 ? 0 : 1   # Offset between vertex-based indices (in `t`)
    sⱼ = δⱼ<0 ? 0 : 1   # and cell-based indices (in `v`)

    vᵢⱼ = v[i-sᵢ, j-sⱼ]
    tᵢ  = t[i-δᵢ, j]
    tⱼ  = t[i, j-δⱼ]

    # Hyp: sign(∇x) == sign(δi) && sign(∇y) == sign(δj)
    #
    # a⋅tᵢⱼ² + b⋅tᵢⱼ + c = 0
    a = 2
    b = -2*(tᵢ + tⱼ)
    c = tᵢ^2 + tⱼ^2 - vᵢⱼ^2
    Δ = b^2 - 4*a*c

    t₁ = t₂ = typemax(T)
    if Δ>0                         # TODO: really need to compute both roots?
        t₁ = (-b + √(Δ)) / 2a
        t₁ > max(tᵢ, tⱼ) || (t₁ = typemax(T))

        t₂ = (-b - √(Δ)) / 2a
        t₂ > max(tᵢ, tⱼ) || (t₂ = typemax(T))
    end

    # Hyp: ∇y == 0  or  ∇x == 0
    t₃ = vᵢⱼ + tᵢ
    t₄ = vᵢⱼ + tⱼ

    # Update if necessary
    min(t₁, t₂, t₃, t₄)
end


const CI = CartesianIndex{2}

struct FastMarching{T}
    t          :: Matrix{T}
    v          :: Matrix{T}
    considered :: PriorityQueue{CI, T, Base.Order.ForwardOrdering}
end

FastMarching(m::Int, n::Int, T=Float64) = FastMarching(
    fill(typemax(T), m+1, n+1),
    fill(0.0, m, n),
    PriorityQueue{CI, T}()
)

function init!(fm::FastMarching, source)
    i = CI(source)
    fm.considered[i] = 0
    fm.t[i] = 0
    fm
end

function march!(fm::FastMarching{T};
                tmax    = typemax(T),
                itermax = typemax(Int),
                verbose = false,
                ) where {T}
    t = fm.t
    v = fm.v
    considered = fm.considered

    @assert size(t) == size(v) .+ (1,1)
    valid = CartesianIndices(t)

    k = 0
    tᵢⱼ = 0.0
    converged = false
    print_progress = 0
    @inbounds while true
        (converged = isempty(considered)) && break
        ij = dequeue!(considered)

        for (δ, dirs) in ((CI( 1,  0), (CI( 1, -1), CI( 1,  1))),
                          (CI( 0,  1), (CI( 1,  1), CI(-1,  1))),
                          (CI(-1,  0), (CI(-1,  1), CI(-1, -1))),
                          (CI( 0, -1), (CI(-1, -1), CI( 1, -1))))
            n = ij + δ
            n ∈ valid || continue
            tₙ = typemax(T)

            for dir in dirs
                (n - dir) ∈ valid || continue
                tₙ = min(tₙ, update(t, v, Tuple(n), Tuple(dir)))
            end

            if tₙ < t[n]
                t[n] = tₙ
                considered[n] = tₙ
            end
        end

        (k += 1) == itermax    && break
        (tᵢⱼ = t[ij]) >= tmax  && break

        if verbose && k > print_progress
            progress = 100 * k / prod(size(t)) |> round |> Int
            print_progress = (progress + 1) * prod(size(t)) / 100
            @printf(stderr, "%3d%%   t=%.2e\r", progress, tᵢⱼ)
        end
    end

    verbose && @printf("%10d %10.2f\n", k, tᵢⱼ)
    converged
end


# Trajectory reconstruction using a gradient descent algorithm
function ray(t::AbstractMatrix{T}, pos; ρ=T(0.1)) where {T}
    # ρ is the step of the gradient descent (in terms of cell size)
    # should be O(1) but can be decreased in cases where the gradient can be very large

    res = [pos]

    (i,j) = pos               # We need to keep track of integer cell indices
    (x,y) = eltype(t).(pos)   # as well as floating-point position

    # Gradient descent
    while true
        tn = get(t, CartesianIndex(i+1, j), 2*t[i,j])
        ts = get(t, CartesianIndex(i-1, j), 2*t[i,j])
        tw = get(t, CartesianIndex(i, j+1), 2*t[i,j])
        te = get(t, CartesianIndex(i, j-1), 2*t[i,j])
        ∇ᵢ = (tn-ts) / 2   # Centered differences
        ∇ⱼ = (tw-te) / 2   # (inconsistent with the FSM: is it a problem?)
        ∇  = (∇ᵢ, ∇ⱼ)

        # t is not differentiable near the origin
        # => stop a bit before arriving
        t[i,j] <= ρ * norm(∇) && break

        # Fixed step length ρ
        x -= ρ * ∇ᵢ / norm(∇); i = round(Int, x)
        y -= ρ * ∇ⱼ / norm(∇); j = round(Int, y)

        # Only record point if it is new
        (i,j) != last(res) && push!(res, (i,j))
    end

    res
end


function from_png(T, filename, colors)
    coeffs = map(colors) do (color, invspeed)
        parse(Colorant, color) => invspeed
    end |> Dict

    img = load(filename)
    sol = FastSweeping(size(img)...)

    col = collect(keys(coeffs))
    map!(sol.v, img) do c
        (_, i) = findmin(c′->colordiff(c,c′), col)
        coeffs[col[i]]
    end

    sol
end

FastSweeping(filename::String, colors) = from_png(FastSweeping, filename, colors)
FastMarching(filename::String, colors) = from_png(FastMarching, filename, colors)


# Convert from vertex-based arrival times to cell-based times
function vertex2cell(t)
    (m,n) = size(t) .- 1

    t′ = similar(t, m, n)
    @inbounds for i in 1:m
        for j in 1:n
            t′[i, j] = min(t[i,j], t[i+1,j], t[i,j+1], t[i+1,j+1])
        end
    end
    t′
end

end
