using IntervalArithmetic

include("./utils.jl")
include("../twofoldInterval/index.jl")
include("../twofoldIntervalPCM/index.jl")

function twofoldIntervalPCM2Weights(
        A::Matrix{TwofoldInterval{T}}
        ):: @NamedTuple{ ŵᴸ::Vector{T}, ŵᵁ::Vector{T} } where {T <: Real}
    m, n = size(A)

    if !isTwofoldIntervalPCM(A)
        throw(ArgumentError("A is not a twofold interval PCM."))
    end

    (Āᴸ, Āᵁ, Āᶜ) = twofoldIntervalPCM2CrispPCM(A)

    Wᶜ = sum(map(j -> geometric_mean(Āᶜ[j,:]), 1:n))
    w̅ᴸ = map(i -> geometric_mean(Āᴸ[i,:]) ./ Wᶜ, 1:n)
    w̅ᵁ = map(i -> geometric_mean(Āᵁ[i,:]) ./ Wᶜ, 1:n)

    ŵᴸ = map(
        i -> inv(max(
            w̅ᴸ[i],
            1 - sum(map(
                j -> w̅ᵁ[j],
                filter(j -> i != j, 1:n)
            )))), 1:n)
    ŵᵁ = map(
        i -> inv(min(
            w̅ᵁ[i],
            1 - sum(map(
                j -> w̅ᴸ[j],
                filter(j -> i != j, 1:n)
            )))), 1:n)

    return (ŵᴸ=ŵᴸ, ŵᵁ=ŵᵁ)
end
