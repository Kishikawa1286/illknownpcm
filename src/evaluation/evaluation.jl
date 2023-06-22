function coincidenceIndex(
    A::Matrix{Interval{T}},
    B::Matrix{Interval{T}}
    )::T where {T <: Real}
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

    conincidenceIndices = fill(1.0, (n, n))
    emptySetMap = fill(false, (n, n))
    for i = 1:n, j = 1:n
        if i == j continue end

        if nearlyEqual(A[i,j].lo, B[i,j].lo) &&
                nearlyEqual(A[i,j].hi, B[i,j].hi)
            conincidenceIndices[i,j] = 1
            continue
        end

        intersection = A[i,j] ∩ B[i,j]
        hull = A[i,j] ∪ B[i,j]

        if !iscommon(intersection) && !nearlyEqual(A[i,j].hi, B[i,j].lo) && !nearlyEqual(A[i,j].lo, B[i,j].hi)
            emptySetMap[i,j] = true
        end

        numerator = iscommon(intersection) ?
            intersection.hi - intersection.lo : 0
        denominator = hull.hi - hull.lo

        # 二つの区間が同じ一点である場合に分母 0 になる
        conincidenceIndices[i,j] = denominator == 0 ? 0 : numerator / denominator
    end

    total = 0.0
    for i = 1:n, j = 1:n
        if i != j && !emptySetMap[i, j]
            total += conincidenceIndices[i, j]
        end
    end

    return total / (n * (n - 1))
end
