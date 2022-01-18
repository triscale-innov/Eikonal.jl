module FastSweeping
using LinearAlgebra
export sweep!, ray


# Fast Sweeping Method
function sweep!(t::AbstractMatrix{T}, v::AbstractMatrix{T} ;
                verbose = false,
                epsilon = eps(T),
                ) where {T}
    m, n = size(v); @assert size(t) == (m+1, n+1)
    change = ones(T, 4)

    ordered_indices(N, δ) = δ>0 ? (2:1:N+1) : (N:-1:1)
    quadrants = ((1,1), (-1,1), (-1,-1), (1, -1))

    # Global iterations
    @inbounds for k in 1:100

        # One sweep per quadrant
        for l in 1:4
            (δᵢ,δⱼ) = quadrants[l]
            maxchange = change[l] = zero(T)
            maximum(change) <= epsilon && return t

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
        end
    end

    error("Max number of iterations reached!")
end


# Trajectory reconstruction using a gradient descent algorithm
function ray(t::AbstractMatrix{T}, pos; ρ=T(0.5)) where {T}
    # ρ is the step of the gradient descent (in terms of cell size)
    # should be O(1) but can be decreased in cases where the gradient can be very large

    res = [pos]

    (i,j) = pos               # We need to keep track of integer cell indices
    (x,y) = eltype(t).(pos)   # as well as floating-point position

    # Gradient descent
    while true
        ∇ᵢ = (t[i+1,j] - t[i-1,j]) / 2   # Centered differences
        ∇ⱼ = (t[i,j+1] - t[i,j-1]) / 2   # (inconsistent with the FSM: is it a problem?)
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

end
