using IntervalArithmetic

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