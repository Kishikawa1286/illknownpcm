# Compute arithmetic mean of non-diagonal elements
function meanOfNonDiagonalElements(A::Matrix{T})::T where {T <: Real}
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
function varianceOfNonDiagonalElements(A::Matrix{T})::T where {T <: Real}
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

function meanOffDiagonalElements(A::Matrix{T})::T where {T <: Real}
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

function P(
    A::Matrix{Interval{T}},
    B::Matrix{Interval{T}}
    )::Matrix{T} where {T <: Real}
    if !isIntervalPCM(A)
        throw(ArgumentError("Given matrix A is not valid as interval PCM."))
    end
    if !isIntervalPCM(B)
        throw(ArgumentError("Given matrix A is not valid as interval PCM."))
    end
    if size(A) != size(B)
        throw(ArgumentError("Given matrices have different size."))
    end

    m, n = size(A)

    indices = fill(1.0, (n, n))
    for i = 1:n, j = 1:n
        if i == j continue end

        if nearlyEqual(A[i,j].lo, B[i,j].lo) &&
                nearlyEqual(A[i,j].hi, B[i,j].hi)
            indices[i,j] = 1
            continue
        end

        intersection = A[i,j] ∩ B[i,j]
        hull = A[i,j] ∪ B[i,j]

        numerator = iscommon(intersection) ?
            intersection.hi - intersection.lo : 0
        denominator = hull.hi - hull.lo

        indices[i,j] = denominator == 0 ? 0 : numerator / denominator
    end

    return indices
end

function Q(
    A::Matrix{Interval{T}},
    B::Matrix{Interval{T}}
    )::Matrix{T} where {T <: Real}
    if !isIntervalPCM(A)
        throw(ArgumentError("Given matrix A is not valid as interval PCM."))
    end
    if !isIntervalPCM(B)
        throw(ArgumentError("Given matrix A is not valid as interval PCM."))
    end
    if size(A) != size(B)
        throw(ArgumentError("Given matrices have different size."))
    end

    m, n = size(A)

    indices = fill(1.0, (n, n))
    emptySetMap = fill(false, (n, n))
    for i = 1:n, j = 1:n
        if i == j continue end

        if !iscommon(B[i,j])
            indices[i,j] = 0
            continue
        end

        if nearlyEqual(A[i,j].lo, B[i,j].lo) &&
                nearlyEqual(A[i,j].hi, B[i,j].hi)
            indices[i,j] = 1
            continue
        end

        intersection = A[i,j] ∩ B[i,j]

        numerator = iscommon(intersection) ?
            intersection.hi - intersection.lo : 0
        denominator = B[i,j].hi - B[i,j].lo

        indices[i,j] = denominator == 0 ? 0 : numerator / denominator
    end

    return indices
end

function R(
    A::Matrix{Interval{T}},
    B::Matrix{Interval{T}}
    )::Matrix{T} where {T <: Real}
    return Q(B, A)
end
