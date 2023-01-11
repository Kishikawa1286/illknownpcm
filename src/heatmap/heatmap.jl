using IntervalArithmetic
using Plots
using LaTeXStrings

include("../nearlyEqual/index.jl")
include("../intervalPCM/index.jl")

function coincidenceIndices(
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

    conincidenceIndices = fill(1.0, (n, n))
    for i = 1:n, j = 1:n
        if i == j continue end

        intersection = A[i,j] ∩ B[i,j]
        hull = A[i,j] ∪ B[i,j]

        numerator = iscommon(intersection) ?
            intersection.hi - intersection.lo : 0
        denominator = hull.hi - hull.lo

        conincidenceIndices[i,j] = numerator / denominator
    end

    return conincidenceIndices
end

function plotConincidenceIndices(
        A::Matrix{Interval{T}},
        B::Matrix{Interval{T}},
        title::LaTeXString
        ) where {T <: Real}
    indices = coincidenceIndices(A, B)
    m, n = size(indices)

    # ヒートマップを [0, 1] スケールにするために表示範囲外に 0 を入れる
    _indices = fill(1.0, (n+1, n+1))
    for i = 1:n+1, j = 1:n+1
        if i != n+1 && j != n+1
            _indices[i,j] = indices[i,j]
        else
            _indices[i,j] = 0.0
        end
    end

    pyplot()

    # 表示範囲外に 0 を入れた行列を使う
    h = heatmap(1:n+1, 1:n+1, _indices,
        c=cgrad([:white, :blue]),
        aspect_ratio=:equal,
        # 表示範囲をヒートマップのタイルに合わせている
        # n+1 は表示しない
        xlims=(0.5,n+0.5), ylims=(0.5,n+0.5),
        xticks=1:n, yticks=1:n,
        # y 軸反転
        yflip=true,
        title=title)
    annotate!(
        [(j, i, text(round(indices[i,j],digits=3),
        12, "DejaVu Serif", :black))
        for i in 1:n for j in 1:n])

    return h
end
