using IntervalArithmetic
using Plots
using LaTeXStrings

include("../nearlyEqual/index.jl")
include("../intervalPCM/index.jl")
include("../nearlyEqual/index.jl")

function coincidenceIndices(
        A::Matrix{Interval{T}},
        B::Matrix{Interval{T}}
        )::Tuple{Matrix{T}, Matrix{Bool}} where {T <: Real}
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

    return (conincidenceIndices, emptySetMap)
end

function plotConincidenceIndices(
        A::Matrix{Interval{T}},
        B::Matrix{Interval{T}},
        title::LaTeXString,
        fontsize::Real=14
        ) where {T <: Real}
    indices, emptySetMap = coincidenceIndices(A, B)
    m, n = size(indices)

    font = Plots.font("sans-serif", fontsize)
    pyplot(guidefont=font, xtickfont=font, ytickfont=font, legendfont=font, colorbar_tickfontsize=fontsize)

    h = heatmap(1:n, 1:n, indices,
        clim=(0, 1),
        c=cgrad([:white, :indianred1]),
        aspect_ratio=:equal,
        # 表示範囲をヒートマップのタイルに合わせている
        xlims=(0.5,n+0.5), ylims=(0.5,n+0.5),
        xticks=1:n, yticks=1:n,
        # y 軸反転
        yflip=true,
        title=title)
    annotate!(
        [(j, i, text(emptySetMap[i,j] ? L"\emptyset" : round(indices[i,j], digits=3),
        fontsize, "sans-serif", :black))
        for i in 1:n for j in 1:n])

    return (h, indices)
end
