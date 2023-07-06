using IntervalArithmetic

# Compute arithmetic mean of non-diagonal elements
@inline function meanOfNonDiagonalElements(A::Matrix{T})::T where {T <: Real}
    m, n = size(A)

    if m != n return 0 end
    if n == 1 return 0 end

    a = 0
    for i = 1:m, j = 1:n
        if i == j continue end
        a += A[i,j]
    end

    return a / n / (n-1)
end

# Compute variance of non-diagonal elements
@inline function varianceOfNonDiagonalElements(A::Matrix{T})::T where {T <: Real}
    m, n = size(A)

    if m != n return 0 end
    if n == 1 return 0 end

    μ = meanOfNonDiagonalElements(A)

    a = 0
    for i = 1:m, j = 1:n
        if i == j continue end
        a += (A[i,j])^2
    end
    a /= n * (n - 1)

    return a - μ^2
end

@inline function comparisonScale(S::Real)
    if S < 1
        throw(ArgumentError("S must be positive int."))
    end
    return vcat(1 ./ (S:-1:2), [1], 2:S)
end

@inline function discretizateIntoComparisonScale(a::T, S::Real)::T where {T <: Real}
    # sorted_scale = [1/9, 1/8, ..., 1/2, 1, 2, ..., 8, 9] if S = 9
    sorted_scale = sort(comparisonScale(S))
    # boundaries = [1/sqrt(9*8), 1/sqrt(8*7), ..., 1/sqrt(3*2), 1/sqrt(2*1), sqrt(1*2), sqrt(2*3), ..., sqrt(7*8), sqrt(8*9)] if S = 9
    boundaries = [sqrt(sorted_scale[i] * sorted_scale[i+1]) for i in 1:length(sorted_scale)-1]

    for i in eachindex(boundaries)
        if a ≤ boundaries[i]
            return sorted_scale[i]
        end
    end

    return last(sorted_scale)
end

@inline function discretizateIntoComparisonScale(A::Matrix{T}, S::Real)::Matrix{T} where {T <: Real}
    return map(a -> discretizateIntoComparisonScale(a, S), A)
end

@inline function discretizateIntoComparisonScale(A::Matrix{Interval{T}}, S::Real)::Matrix{Interval{T}} where {T <: Real}
    return map(a -> discretizateIntoComparisonScale(a.lo, S)..discretizateIntoComparisonScale(a.hi, S), A)
end
