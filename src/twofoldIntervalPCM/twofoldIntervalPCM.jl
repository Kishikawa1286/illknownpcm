using IntervalArithmetic

include("../intervalPCM/index.jl")
include("../twofoldInterval/index.jl")

@inline function isTwofoldIntervalPCM(
        A::Matrix{TwofoldInterval{T}})::Bool where {T <: Real}
    m, n = size(A)
    # PCMは正方行列
    if m != n return false end

    for i = 1:n, j = 1:n
        # 外側Aᵢⱼ⁺が空集合の成分が存在しないか
        if !iscommon(A[i,j][2]) return false end

        # 内側Aᵢⱼ⁻が空集合の場合
        if !iscommon(A[i,j][1])
            # 対角成分は([1,1], [1,1])
            if i == j return false end

            # (i,j)成分の外の両端
            aᵢⱼᴸ⁺ = A[i,j][2].lo; aᵢⱼᵁ⁺ = A[i,j][2].hi
            # (j,i)成分の外の両端
            aⱼᵢᴸ⁺ = A[j,i][2].lo; aⱼᵢᵁ⁺ = A[j,i][2].hi

            # 外側のreciprocity
            if !nearlyEqual(aᵢⱼᴸ⁺, 1 / aⱼᵢᵁ⁺) return false end
            if !nearlyEqual(aᵢⱼᵁ⁺, 1 / aⱼᵢᴸ⁺) return false end

            continue
        end

        # (i,j)成分の内の両端
        aᵢⱼᴸ⁻ = A[i,j][1].lo; aᵢⱼᵁ⁻ = A[i,j][1].hi
        # (i,j)成分の外の両端
        aᵢⱼᴸ⁺ = A[i,j][2].lo; aᵢⱼᵁ⁺ = A[i,j][2].hi

        # aᵢⱼᴸ⁺ ≤ aᵢⱼᴸ⁻ ≤ aᵢⱼᵁ⁻ ≤ aᵢⱼᵁ⁺ を満たすか
        if aᵢⱼᴸ⁺ > aᵢⱼᴸ⁻ return false end
        if aᵢⱼᴸ⁻ > aᵢⱼᵁ⁻ return false end
        if aᵢⱼᵁ⁻ > aᵢⱼᵁ⁺ return false end

        # 対角成分が([1,1], [1,1])であるか
        # ≈ 1 は許容
        if i == j
            if !nearlyEqual(aᵢⱼᴸ⁻, 1) return false end
            if !nearlyEqual(aᵢⱼᵁ⁻, 1) return false end
            if !nearlyEqual(aᵢⱼᴸ⁺, 1) return false end
            if !nearlyEqual(aᵢⱼᵁ⁺, 1) return false end
        end

        # (j,i)成分の内の両端
        aⱼᵢᴸ⁻ = A[j,i][1].lo; aⱼᵢᵁ⁻ = A[j,i][1].hi
        # (j,i)成分の外の両端
        aⱼᵢᴸ⁺ = A[j,i][2].lo; aⱼᵢᵁ⁺ = A[j,i][2].hi

        # reciprocityが成り立つか
        if !nearlyEqual(aᵢⱼᴸ⁻, 1 / aⱼᵢᵁ⁻) return false end
        if !nearlyEqual(aᵢⱼᵁ⁻, 1 / aⱼᵢᴸ⁻) return false end
        if !nearlyEqual(aᵢⱼᴸ⁺, 1 / aⱼᵢᵁ⁺) return false end
        if !nearlyEqual(aᵢⱼᵁ⁺, 1 / aⱼᵢᴸ⁺) return false end
    end

    return true
end

"""
Interval PCM 2 個から Twofold Interval PCM を生成する
"""
function intervalPCM2TwofoldIntervalPCM(
        A₁::Matrix{Interval{T}},
        A₂::Matrix{Interval{T}}
        )::Matrix{TwofoldInterval{T}} where {T <: Real}
    if !isIntervalPCM(A₁)
        throw(ArgumentError("A₁ is not valid as interval PCM."))
    end
    if !isIntervalPCM(A₂)
        throw(ArgumentError("A₂ is not valid as interval PCM."))
    end

    if size(A₁) != size(A₂)
        throw(ArgumentError("Size of A₁ and A₂ are differnt."))
    end

    m, n = size(A₁)
    # ([1, 1], [1, 1])で埋める
    A = fill((1..1, 1..1), (n, n))

    for i = 1:n, j = 1:n
        # 対角成分は初期値の([1, 1], [1, 1])のままで良いのでスキップ
        if i == j continue end

        # Aᵢⱼ⁻はA₁ᵢⱼとA₂ᵢⱼの共通部分
        # Aᵢⱼ⁺はA₁ᵢⱼとA₂ᵢⱼの和集合の凸包
        #　∪演算子は和集合でなく和集合の凸包を返す
        A[i,j] = (A₁[i,j] ∩ A₂[i,j], A₁[i,j] ∪ A₂[i,j])
    end

    if !isTwofoldIntervalPCM(A)
        throw(ErrorException("Calculation Result A is not a PCM."))
    end

    return A
end
