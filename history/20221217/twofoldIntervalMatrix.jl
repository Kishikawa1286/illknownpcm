using LinearAlgebra
using IntervalArithmetic
using Printf

include("utils.jl")

# (Aᵢⱼ⁻, Aᵢⱼ⁺)
TwofoldInterval = Tuple{Interval{T}, Interval{T}} where {T <: Real}

@inline function isTwofoldIntervalPCM(
        A::Matrix{TwofoldInterval{T}}
        )::Bool where {T <: Real}
    # 正方行列でない
    m, n = size(A)
    if m != n
        return false
    end

    for i = 1:n, j = 1:n
        aᵢⱼᴸ⁺ = A[i,j][2].lo
        aᵢⱼᴸ⁻ = A[i,j][1].lo
        aᵢⱼᵁ⁻ = A[i,j][1].hi
        aᵢⱼᵁ⁺ = A[i,j][2].hi
        # aᵢⱼᴸ⁺ ≤ aᵢⱼᴸ⁻ ≤ aᵢⱼᵁ⁻ ≤ aᵢⱼᵁ⁺ を満たさない
        if !(aᵢⱼᴸ⁺ <= aᵢⱼᴸ⁻ && aᵢⱼᴸ⁻ <= aᵢⱼᵁ⁻ && aᵢⱼᵁ⁻ <= aᵢⱼᵁ⁺)
            return false
        end
        # 対角成分が aᵢⱼᴸ⁺ = aᵢⱼᴸ⁻ = aᵢⱼᵁ⁻ = aᵢⱼᵁ⁺ = 1 を満たさない
        if i == j && (aᵢⱼᴸ⁺ != 1 || aᵢⱼᵁ⁺ != 1 || aᵢⱼᴸ⁻ != 1 || aᵢⱼᵁ⁻ != 1)
            return false
        end
        aⱼᵢᴸ⁺ = A[j,i][2].lo
        aⱼᵢᴸ⁻ = A[j,i][1].lo
        aⱼᵢᵁ⁻ = A[j,i][1].hi
        aⱼᵢᵁ⁺ = A[j,i][2].hi
        # reciprocity aᵢⱼᴸ⁺ = 1/aⱼᵢᵁ⁺, aᵢⱼᴸ⁻ = 1/aⱼᵢᵁ⁻ を満たさない
        if !(
                aᵢⱼᴸ⁺ * aⱼᵢᵁ⁺ - 1 < 1e-5 &&
                aᵢⱼᴸ⁻ * aⱼᵢᵁ⁻ - 1 < 1e-5 &&
                aⱼᵢᴸ⁺ * aᵢⱼᵁ⁺ - 1 < 1e-5 &&
                aⱼᵢᴸ⁻ * aᵢⱼᵁ⁻ - 1 < 1e-5)
            return false
        end
    end

    return true
end

function generateTwofoldIntervalMatrix(
        A₁::Matrix{Interval{T}},
        A₂::Matrix{Interval{T}}
        )::Matrix{TwofoldInterval{T}} where {T <: Real}
    if size(A₁) != size(A₂)
        throw(ArgumentError("""
                Size of A₁ and A₂ are differnt.
                """))
    end

    m, n = size(A₁) # 先頭の行列の大きさ（すべての行列の大きさは同じ）
    A = [(1..1, 1..1) for a in 1:n, b in 1:n] # m*nを(Aᵢⱼ⁺ = [1, 1], Aᵢⱼ⁻ = [1, 1])で埋める

    for i = 1:n, j = 1:n
        if i == j
            continue # 対角成分は初期値のままで良いのでスキップ
        end

        # Aᵢⱼ⁻はA₁ᵢⱼとA₂ᵢⱼの共通部分
        # Aᵢⱼ⁺はA₁ᵢⱼとA₂ᵢⱼの和集合の凸包
        #　∪演算子は和集合でなく和集合の凸包を返す
        A[i,j] = (A₁[i,j] ∩ A₂[i,j], A₁[i,j] ∪ A₂[i,j])
    end

    return A
end

# A ⊆ B なら true
function allIncluded(A::Matrix{TwofoldInterval{T}},
        B::Matrix{TwofoldInterval{T}}, k::Number) where {T <: Real}
    if size(A) != size(B)
        return false
    end
    m, n = size(A)
    for i = 1:m, j = 1:n
        if iscommon(A[i,j][k]) && iscommon(B[i,j][k])
            if !(A[i,j][k].lo >= B[i,j][k].lo && A[i,j][k].hi <= B[i,j][k].hi)
                return false
            end
        else
            return false
        end
        if iscommon(A[i,j][k]) && !iscommon(B[i,j][k])
            return false
        end
    end
    return true
end

function nearlyEqualMatrix(A::Matrix{TwofoldInterval{T}},
        B::Matrix{TwofoldInterval{T}}) where {T <: Real}
    if size(A) != size(B)
        return false
    end
    m, n = size(A)
    for i = 1:m, j = 1:n
        # 片方だけi,j成分が空集合
        if iscommon(A[i,j][1]) != iscommon(B[i,j][1])
            return false
        end
        # 内の空集合以外
        if iscommon(A[i,j][1]) && iscommon(B[i,j][1])
            if !nearlyEqual(A[i,j][1].lo, B[i,j][1].lo)
                return false
            end
            if !nearlyEqual(A[i,j][1].hi, B[i,j][1].hi)
                return false
            end
        end
        # 外
        if !nearlyEqual(A[i,j][2].lo, B[i,j][2].lo)
            return false
        end
        if !nearlyEqual(A[i,j][2].hi, B[i,j][2].hi)
            return false
        end
    end
    return true
end

function nearlyEqualMatrix(A::Matrix{TwofoldInterval{T}},
        B::Matrix{TwofoldInterval{T}}, k::Number) where {T <: Real}
    if size(A) != size(B)
        return false
    end
    m, n = size(A)
    for i = 1:m, j = 1:n
        # 片方だけi,j成分が空集合
        if iscommon(A[i,j][k]) != iscommon(B[i,j][k])
            return false
        end
        # 内の空集合以外
        if iscommon(A[i,j][k]) && iscommon(B[i,j][k])
            if !nearlyEqual(A[i,j][k].lo, B[i,j][k].lo)
                return false
            end
            if !nearlyEqual(A[i,j][k].hi, B[i,j][k].hi)
                return false
            end
        else
            return false
        end
    end
    return true
end

# ([aᵢⱼᴸ⁻, aᵢⱼᵁ⁻], [aᵢⱼᴸ⁺, aᵢⱼᵁ⁺])
function showElements(A::Matrix{TwofoldInterval{T}}) where {T <: Real}
    m, n = size(A)
    for i = 1:m, j = 1:n
        if iscommon(A[i,j][1])
            print(@sprintf("([%0.3f, %0.3f], [%0.3f, %0.3f])", A[i,j][1].lo, A[i,j][1].hi, A[i,j][2].lo, A[i,j][2].hi))
        else
            print(@sprintf("(      ∅       , [%0.3f, %0.3f])", A[i,j][2].lo, A[i,j][2].hi))
        end
        if j == n
            print('\n')
        else
            print(" ")
        end
    end
end

# [ aᵢⱼᴸ⁺, [aᵢⱼᴸ⁻, aᵢⱼᵁ⁻], aᵢⱼᵁ⁺]
function showElements2(A::Matrix{TwofoldInterval{T}}) where {T <: Real}
    m, n = size(A)
    for i = 1:m, j = 1:n
        if iscommon(A[i,j][1])
            print(@sprintf("[%0.3f, [%0.3f, %0.3f], %0.3f]", A[i,j][2].lo, A[i,j][1].lo, A[i,j][1].hi, A[i,j][2].hi))
        else
            print(@sprintf("[%0.3f,       ∅       , %0.3f]", A[i,j][2].lo, A[i,j][2].hi))
        end
        if j == n
            print('\n')
        else
            print(" ")
        end
    end
end

function showWidths(A::Matrix{TwofoldInterval{T}}) where {T <: Real}
    m, n = size(A)
    for i = 1:m, j = 1:n
        if iscommon(A[i,j][1])
            print(@sprintf("(%0.5f, %0.5f)", A[i,j][1].hi - A[i,j][1].lo, A[i,j][2].hi - A[i,j][2].lo))
        else
            print(@sprintf("( undef , %0.5f)", A[i,j][2].hi - A[i,j][2].lo))
        end
        if j == n
            print('\n')
        else
            print(" ")
        end
    end
end
