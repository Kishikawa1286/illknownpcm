using IntervalArithmetic

include("./method3.jl")
include("../nearlyEqual/index.jl")
include("../twofoldInterval/index.jl")
include("../twofoldIntervalPCM/index.jl")

function generateRevisedMatrix(
        Aₖ::Matrix{Interval{T}}
        )::Tuple{Matrix{Interval{T}}, Matrix{Interval{T}}} where {T <: Real}
    lpResult = solveApproximationLP_m3(Aₖ)
    wₖᴸ⁻ = lpResult.wₖᴸ⁻; wₖᵁ⁻ = lpResult.wₖᵁ⁻
    wₖᴸ⁺ = lpResult.wₖᴸ⁺; wₖᵁ⁺ = lpResult.wₖᵁ⁺
    n = length(wₖᴸ⁻)

    Ãₖ⁻ = fill(1.0..1.0, (n, n)); Ãₖ⁺ = fill(1.0..1.0, (n, n))
    for i = 1:n, j = 1:n
        if i == j continue end

        ãₖᵢⱼᴸ⁺ = wₖᴸ⁺[i] / wₖᵁ⁺[j]
        ãₖᵢⱼᵁ⁺ = wₖᵁ⁺[i] / wₖᴸ⁺[j]
        Ãₖ⁺[i,j] = ãₖᵢⱼᴸ⁺..ãₖᵢⱼᵁ⁺

        if any(isinf.([wₖᴸ⁻[i], wₖᵁ⁻[i], wₖᵁ⁻[j], wₖᴸ⁻[j]]))
            Ãₖ⁻[i,j] = emptyinterval()
            continue
        end

        ãₖᵢⱼᴸ⁻ = wₖᴸ⁻[i] / wₖᵁ⁻[j]
        ãₖᵢⱼᵁ⁻ = wₖᵁ⁻[i] / wₖᴸ⁻[j]
        Ãₖ⁻[i,j] = ãₖᵢⱼᴸ⁻..ãₖᵢⱼᵁ⁻
    end

    return (Ãₖ⁻, Ãₖ⁺)
end

