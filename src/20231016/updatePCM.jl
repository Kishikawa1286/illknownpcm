using IntervalArithmetic

include("./solution.jl")
include("../twofoldInterval/index.jl")
include("../twofoldIntervalPCM/index.jl")

function updatePCM(
        A::Matrix{TwofoldInterval{T}},
        solution::Solution{T}
        )::Matrix{TwofoldInterval{T}} where {T <: Real}
    if !isTwofoldIntervalPCM(A)
        throw(ArgumentError("Given matrix is not valid as twofold interval matrix."))
    end

    m, n = size(A)
    Â = deepcopy(A)

    wᴸ = solution.wᴸ; wᵁ = solution.wᵁ
    wᴸ⁻ = solution.wᴸ⁻; wᵁ⁻ = solution.wᵁ⁻
    wᴸ⁺ = solution.wᴸ⁺; wᵁ⁺ = solution.wᵁ⁺

    for i = 1:n, j = 1:n
        # 対角成分は (1..1, 1..1) で固定なので更新不要
        if i == j continue end

        aᵢⱼᴸ⁻ = A[i,j][1].lo; aᵢⱼᵁ⁻ = A[i,j][1].hi
        aᵢⱼᴸ⁺ = A[i,j][2].lo; aᵢⱼᵁ⁺ = A[i,j][2].hi
        wᵢᴸ = wᴸ[i]; wᵢᵁ = wᵁ[i]
        wᵢᴸ⁻ = wᴸ⁻[i]; wᵢᵁ⁻ = wᵁ⁻[i]; wᵢᴸ⁺ = wᴸ⁺[i]; wᵢᵁ⁺ = wᵁ⁺[i]
        wⱼᴸ = wᴸ[j]; wⱼᵁ = wᵁ[j]
        wⱼᴸ⁻ = wᴸ⁻[j]; wⱼᵁ⁻ = wᵁ⁻[j]; wⱼᴸ⁺ = wᴸ⁺[j]; wⱼᵁ⁺ = wᵁ⁺[j]

        âᵢⱼᴸ⁺ = min(aᵢⱼᴸ⁺, wᵢᴸ⁻/wⱼᵁ, wᵢᴸ/wⱼᵁ⁻)
        âᵢⱼᴸ⁻ = max(aᵢⱼᴸ⁻, wᵢᴸ⁺/wⱼᵁ, wᵢᴸ/wⱼᵁ⁺)
        âᵢⱼᵁ⁻ = min(aᵢⱼᵁ⁻, wᵢᵁ⁺/wⱼᴸ, wᵢᵁ/wⱼᴸ⁺)
        âᵢⱼᵁ⁺ = max(aᵢⱼᵁ⁺, wᵢᵁ⁻/wⱼᴸ, wᵢᵁ/wⱼᴸ⁻)

        # âᵢⱼᵁ⁺ = âᵢⱼᴸ⁺ の場合などに precision error で âᵢⱼᵁ⁺ < âᵢⱼᴸ⁺ となることがある
        # âᵢⱼᴸ⁺ と âᵢⱼᵁ⁺ が十分に近い値ならば âᵢⱼᴸ⁺ <- âᵢⱼᵁ⁺
        âᵢⱼᴸ⁺ = min(âᵢⱼᴸ⁺, âᵢⱼᵁ⁺)
        # âᵢⱼᴸ⁻ と âᵢⱼᵁ⁻ が十分に近い値ならば âᵢⱼᴸ⁻ <- âᵢⱼᵁ⁻
        âᵢⱼᴸ⁻ = min(âᵢⱼᴸ⁻, âᵢⱼᵁ⁻)

        # (Âᵢⱼ⁻, Âᵢⱼ⁺)
        if âᵢⱼᴸ⁻ > âᵢⱼᵁ⁻
            Â[i, j] = (emptyinterval(), âᵢⱼᴸ⁺..âᵢⱼᵁ⁺)
        else
            Â[i, j] = (âᵢⱼᴸ⁻..âᵢⱼᵁ⁻, âᵢⱼᴸ⁺..âᵢⱼᵁ⁺)
        end
    end

    return Â
end
