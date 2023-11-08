Solution = @NamedTuple{
    wᴸ::Vector{T}, wᵁ::Vector{T},
    wᴸ⁻::Vector{T}, wᵁ⁻::Vector{T},
    wᴸ⁺::Vector{T}, wᵁ⁺::Vector{T},
    εᴸ::Vector{T}, εᵁ::Vector{T},
    optimalValue::T
    } where {T <: Real}
