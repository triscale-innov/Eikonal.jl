module Eikonal

using DataStructures
using Images
using LinearAlgebra
using Printf

export FastSweeping, sweep!
export FastMarching, init!, march!
export ray
export vertex2cell


"""
    brgc(n)

Get the Binary-Reflected Gray Code list for `n` bits (most significant bit last,
i.e. using the reflect-and-suffix method).
"""
brgc(n::Integer) = brgc(Val(n))
brgc(::Val{1}) = [(0,), (1,)]
function brgc(::Val{N}) where {N}
    r = brgc(Val(N-1))
    res1 = map(r) do c
        (c..., 0)
    end
    res2 = map(reverse(r)) do c
        (c..., 1)
    end
    return append!(res1, res2)
end

struct Orthant{D}
    δ :: NTuple{D, Int}
    s :: NTuple{D, Int}
    rev :: NTuple{D, Bool}

    function Orthant(code::NTuple{D, Int}) where {D}
        δ = map(code) do cᵢ
            cᵢ == 0 ? 1 : -1
        end
        s = map(δ) do δᵢ
            δᵢ < 0 ? 0 : 1
        end
        rev = map(δ) do δᵢ
            δᵢ < 0 ? true : false
        end
        new{D}(δ, s, rev)
    end
end

struct FastSweeping{T, D}
    t      :: Array{T, D}
    v      :: Array{T, D}
    change :: Vector{T}
    iter   :: Base.RefValue{Tuple{Int, Int}}  # (k, l)

    function FastSweeping(σ::Array{T, D}) where {T,D}
        siz = size(σ)
        t = fill(typemax(T), (siz .+ 1)...)
        change = ones(T, 2^D)
        iter = Ref((1, 0))
        new{T,D}(t, σ, change, iter)
    end
end

function FastSweeping(T::DataType, siz::NTuple{D, Int}) where {D}
    σ = fill(typemax(T), siz...)
    FastSweeping(σ)
end

function Base.show(io::IO, fs::FastSweeping)
    grid_size = join(string.(size(fs.v)), "×")
    print(io, "FastSweeping solver on a $grid_size grid")
end

function init!(fs::FastSweeping, source)
    fs.t[source...] = 0
    fs.change .= Inf
    fs
end

"""
    subtuples(t::NTuple{N, T}) where {N, T}

Generate the list of all sub-tuples obtained by removing one element from `t`.
"""
@generated function subtuples(t::NTuple{N, T}) where {N, T}
    expr = Expr(:tuple)
    expr.args = map(1:N) do i
        Expr(:tuple, (:(t[$j]) for j in 1:N if i != j)...)
    end
    expr
end

@inline function update(tᵢ::NTuple{N, T}, v) where {N, T}
    square(x) = x*x

    # Hyp: ∇t is in an N-dimensional orthant
    #      tᵢ contains the values of t in the N neighboring vertices
    #
    # t satisfies:
    #    a⋅t² + b⋅t + c = 0
    a = N
    b = -2*sum(tᵢ)
    c = sum(square, tᵢ) - v^2
    Δ = b^2 - 4*a*c

    t = typemax(T)
    if Δ >= 0
        t = (-b + √(Δ)) / 2a
        t <= maximum(tᵢ) && (t = typemax(T))
    end

    # Hyp: ∇t lies on the boundary of the orthant
    # => perform N-1 dimensional updates
    for tᵢ′ in subtuples(tᵢ)
        t = min(t, update(tᵢ′, v))
    end
    return t
end

@inline function update(tᵢ::NTuple{1, T}, v) where {T}
    return tᵢ[1] + v
end

# Fast Sweeping Method
function sweep!(fs :: FastSweeping{T};
                verbose = false,
                epsilon = eps(T),
                nsweeps = typemax(Int),
                ) where {T}
    D = ndims(fs.v)

    square(x) = x*x
    function indices(axis, rev)
        a = rev ? reverse(axis) : StepRange(axis)
        a[2:end]
    end

    # Build a list of all orthants, in the order determined by a Gray code
    # (since it forms a Hamiltonian cycle on the N-D hypercube)
    orthants = Orthant.(brgc(D))

    k = first(fs.iter[]) # Global iter number
    l = last(fs.iter[])  # Quadrant index
    sweep = 0            # Sweep number
    @inbounds while true
        l += 1
        sweep += 1
        if l > length(orthants)
            l = 1
            k += 1
        end
        fs.iter[] = (k, l)

        maxchange = fs.change[l] = zero(T)
        maximum(fs.change) <= epsilon && return true

        orthant = orthants[l]
        δ = orthant.δ     # direction (±1) alongside each axis
        s = orthant.s     # shift/offset between vertex-based indices (in `t`)
                          #    and cell-based indices (in `v`)
        rev = orthant.rev # rev: indicates whether each axis should be reversed

        for I in CartesianIndices(indices.(axes(fs.t), rev))
            v = fs.v[I - CartesianIndex(s)]

            tᵢ = ntuple(Val(D)) do k
                dir = ntuple(l -> k==l ? δ[k] : 0, Val(D))
                fs.t[I - CartesianIndex(dir)]
            end
            t = update(tᵢ, v)

            # Update if necessary
            if t < fs.t[I]
                maxchange = max(maxchange, (fs.t[I]-t)/t)
                fs.t[I] = t
            end
        end
        verbose && println("iter $k, sweep $l: change = $maxchange")
        fs.change[l] = maxchange

        sweep >= nsweeps && return false
    end
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
ray(t::AbstractArray{T,D}, pos) where {T,D} = ray(t, pos, Integration(0.1))

struct Integration
    rho :: Float64
end

function ray(t::AbstractArray{T,D}, pos, method::Integration) where {T,D}
    # ρ is the step of the gradient descent (in terms of cell size)
    # should be O(1) but can be decreased in cases where the gradient can be very large
    ρ = method.rho

    res = [pos]

    # We need to keep track of floating point coordinates
    x = eltype(t).(pos)

    # Gradient descent
    while true
        I = CartesianIndex(round.(Int, x)) # Integer coordinates
        ∇ = ntuple(Val(D)) do k
            δ = CartesianIndex(ntuple(l -> k==l ? 1 : 0, Val(D)))
            t1 = min(t[I + δ], 2*t[I])
            t2 = min(t[I - δ], 2*t[I])
            (t1-t2) / 2
        end

        t_min = ntuple(Val(D)) do k
            δ = CartesianIndex(ntuple(l -> k==l ? 1 : 0, Val(D)))
            min(t[I + δ], t[I - δ])
        end

        # Stop when we reached a local minimum
        minimum(t_min) >= t[I] && break

        # Update position with a fixed step length ρ
        x = x .- ρ .* ∇ ./ norm(∇)

        # Only record point if it is new
        pos = Tuple(I)
        pos != last(res) && push!(res, pos)
    end

    res
end

struct NearestMin end
function ray(t::AbstractArray{T,D}, pos, ::NearestMin) where {T, D}
    I = CartesianIndex(pos)
    res = [pos]

    while true
        # Find minimal time in a box surrounding the current position
        box = I-oneunit(I):I+oneunit(I)
        _, ibox = findmin(box) do I′
            I′ ∈ CartesianIndices(t) ? t[I′] : Inf
        end
        I′ = box[ibox]

        # Stop when we reached a local minimum
        I == I′ && break

        # Update current position
        I = I′
        push!(res, Tuple(I))
    end

    res
end


function from_png(T, filename, colors)
    σ = from_png(filename, colors)
    T(σ)
end

function from_png(filename::AbstractString, colors)
    img = load(filename)
    dict = map(colors) do (color, invspeed)
        parse(Colorant, color) => invspeed
    end |> Dict

    col = collect(keys(dict))
    map(img) do c
        (_, i) = findmin(c′->colordiff(c,c′), col)
        dict[col[i]]
    end
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

include("precompile.jl")

end
