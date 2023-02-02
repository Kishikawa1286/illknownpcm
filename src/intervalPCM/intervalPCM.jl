using IntervalArithmetic
using Random
using Distributions

include("../nearlyEqual/index.jl")

@inline function isIntervalPCM(A::Matrix{Interval{T}})::Bool where {T <: Real}
    m, n = size(A)
    # PCMは正方行列
    if m != n return false end

    for i = 1:n, j = 1:n
        if !iscommon(A[i,j])
            if iscommon(A[j,i]) return false end
            continue
        end
        
        # (i,j)成分の両端
        aᵢⱼᴸ = A[i,j].lo; aᵢⱼᵁ = A[i,j].hi

        # 対角成分が[1,1]であるか
        if i == j
            if !nearlyEqual(aᵢⱼᴸ, 1) return false end
            if !nearlyEqual(aᵢⱼᵁ, 1) return false end
        end

        # (j,i)成分の両端
        aⱼᵢᴸ = A[j,i].lo; aⱼᵢᵁ = A[j,i].hi

        # reciprocityが成り立つか
        if !nearlyEqual(aᵢⱼᴸ, 1 / aⱼᵢᵁ) return false end
        if !nearlyEqual(aᵢⱼᵁ, 1 / aⱼᵢᴸ) return false end
    end

    return true 
end

function isWeaklyConsistentIntervalPCM(
        A::Matrix{T})::Bool where {T <: Real}
    if !isIntervalPCM(A) return false end

    _, n = size(A)

    for i = 1:n, j = 1:n, k = 1:n, s = 1:n, r = 1:n
        if j == i || k == i || s == i || r == i
            continue
        end
        if j == k || s == r
            continue
        end

        if !nearlyEqualLoose(A[i,j].lo * A[j,k].hi * A[k,i].lo,
                A[i,s].lo * A[s,r].hi * A[r,i].lo)
            return false
        end

        if !nearlyEqual(A[i,j].lo * A[j,k].hi * A[k,i].lo, 1) &&
                A[i,j].lo * A[j,k].hi * A[k,i].lo > 1
            return false
        end

        if !nearlyEqual(A[s,i].lo, A[s,r].hi * A[r,i].lo) &&
                A[s,i].lo > A[s,r].hi * A[r,i].lo
            return false
        end
    end

    return true
end

function isConsistentIntervalPCM(
        A::Matrix{T})::Bool where {T <: Real}
    if !isIntervalPCM(A) return false end
    if !isWeaklyConsistentIntervalPCM(A) return false end

    _, n = size(A)

    for i = 1:n, j = 1:n, k = 1:n
        if j == i || k == i || j == k continue end

        α = A[i,j].lo * A[j,k].hi * A[k,i].lo +
            A[j,k].hi * A[k,i].lo +
            sum(map(r -> A[r,j].hi * A[j,i].lo, filter(r -> r != i && r != j, 1:n)))
        β = minimum(map(
            r -> A[r,i].lo + 1 +
            sum(map(s -> A[s,r].hi * A[r,i].lo, filter(s -> s != i && s != r, 1:n))),
            filter(r -> r != i, 1:n)))
        γ = 1 + sum(map(r -> A[r,i].lo, filter(r -> r != i, 1:n)))
        δ = A[j,k].hi * A[k,i].lo + A[i,j].lo * A[j,k].hi * A[k,i].lo +
            sum(map(r -> A[r,i].lo, filter(r -> r != i && r != j, 1:n)))
        ϵ = maximum(map(r -> A[r,j].hi * A[j,i].lo + A[i,j].lo * A[j,k].hi * A[k,i].lo +
            sum(map(s -> A[s,i].lo, filter(s -> s != i && s != r, 1:n))),
            filter(r -> r != i && r != j, 1:n)))

        if min(α, β, γ) < max(δ, ϵ)
            return false
        end
    end

    return true
end

@inline function randamizedIntervalPCM(
        A::Matrix{T},
        seed::Integer=10,
        width::Real=0.05 # 自然対数スケールでの幅
        )::Matrix{Interval{T}} where {T <: Real}
    # seed 固定
    Random.seed!(seed)

    m, n = size(A)

    if m != n
        throw(ArgumentError("A must be square matrix."))
    end

    B = log.(A)
    C = fill(1..1, (m, n))

    for i = 1:n
        for j = i:n
            if i == j continue end
            r₁ = rand(Uniform(0, width))
            r₂ = rand(Uniform(0, width))
            C[i,j] = exp(B[i,j] - r₁)..exp(B[i,j] + r₂)
        end
    end

    for i = 1:n
        for j = 1:i
            if i == j continue end
            C[i,j] = 1/C[j,i]
        end
    end

    return C
end
