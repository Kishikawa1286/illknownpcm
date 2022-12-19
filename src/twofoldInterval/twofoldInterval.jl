using IntervalArithmetic

include("../nearlyEqual/nearlyEqual.jl")

TwofoldInterval = Tuple{Interval{T}, Interval{T}} where {T <: Real}

@inline function isTwofoldInterval(
        A::TwofoldInterval{T}
        )::Bool where {T <: Real}
    # 外が空集合は不可
    if !iscommon(A[2]) return false end

    # この時点で内が∅ならTwofoldInterval
    if !iscommon(A[1]) return true end

    aᴸ⁻ = A[1].lo; aᵁ⁻ = A[1].hi
    aᴸ⁺ = A[2].lo; aᵁ⁺ = A[2].hi

    # aᴸ⁻ ≥ aᴸ⁺ かどうか検証
    # aᴸ⁻ ≈ aᴸ⁺ は許容する
    if aᴸ⁻ < aᴸ⁺ && !nearlyEqual(aᴸ⁻, aᴸ⁺) return false end

    # aᵁ⁻ ≤ aᵁ⁺ かどうか検証
    # aᵁ⁻ ≈ aᵁ⁺ は許容する
    if aᵁ⁻ > aᵁ⁺ && !nearlyEqual(aᵁ⁻, aᵁ⁺) return false end

    return true
end
