module FastSweeping
using LinearAlgebra
export sweep!, ray


# Fast Sweeping Method
function sweep!(t, v)
    T = eltype(t)
    m, n = size(t)
    change = ones(T, 4)

    ordered_indices(N, δ) = δ>0 ? (2:N) : (N-1:-1:1)
    quadrants = ((1,1), (-1,1), (-1,-1), (1, -1))

    # Global iterations
    @inbounds for k in 1:100

        # One sweep per quadrant
        for l in 1:4
            (δᵢ,δⱼ) = quadrants[l]
            maxchange = change[l] = zero(T)
            maximum(change) <= eps(T) && return t  # TODO: really need to wait for 3 useless sweeps?

            for j in ordered_indices(n, δⱼ)
                for i in ordered_indices(m, δᵢ)
                    vᵢⱼ = v[i,j]                   # TODO: fix off-by-one error
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
            println("iter $k, sweep $l: change = $maxchange")
            change[l] = maxchange
        end
    end

    error("Max number of iterations reached!")
end


# Trajectory reconstruction using a gradient descent algorithm
function ray(t, pos)
    res = [pos]

    # Step (in terms of cell size)
    # should be O(1)
    ρ = eltype(t)(0.5)

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
