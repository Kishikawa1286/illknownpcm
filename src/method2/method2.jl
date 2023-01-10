using JuMP
import HiGHS

include("../nearlyEqual/index.jl")
include("../twofoldInterval/index.jl")
include("../twofoldIntervalPCM/index.jl")

ApproximationLPResult_m2 = @NamedTuple{
    wᴸ⁻::Vector{T}, wᵁ⁻::Vector{T},
    wᴸ⁺::Vector{T}, wᵁ⁺::Vector{T},
    εᴸ::Vector{T}, εᵁ::Vector{T},
    vᴸ⁻::Vector{T}, vᵁ⁻::Vector{T},
    vᴸ⁺::Vector{T}, vᵁ⁺::Vector{T},
    optimalValue::T
    } where {T <: Real}

function solveApproximationLP_m2(
        pcms::Vector{Matrix{Interval{T}}}
        )::ApproximationLPResult_m2{T} where {T <: Real}
    ε = 1e-8

    for k = eachindex(pcms)
        Aₖ = pcms[k]
        if !isIntervalPCM(Aₖ)
            throw(ArgumentError("Matrix k = $(k) is not valid as PCM."))
        end
    end

    m = length(pcms)
    _, n = size(pcms[1])
    model = Model(HiGHS.Optimizer)
    set_silent(model)

    try
        # wᵢᴸ⁻ ≥ ε, wᵢᵁ⁻ ≥ ε
        @variable(model, wᴸ⁻[i=1:n] ≥ ε); @variable(model, wᵁ⁻[i=1:n] ≥ ε)
        # wᵢᴸ⁺ ≥ ε, wᵢᵁ⁺ ≥ ε
        @variable(model, wᴸ⁺[i=1:n] ≥ ε); @variable(model, wᵁ⁺[i=1:n] ≥ ε)
        # εᵢᴸ ≥ 0, εᵢᵁ ≥ 0
        @variable(model, εᴸ[i=1:n] ≥ 0); @variable(model, εᵁ[i=1:n] ≥ 0)
        # vᵢᴸ⁻ ≥ 0, vᵢᵁ⁻ ≥ 0
        @variable(model, vᴸ⁻[i=1:n] ≥ 0); @variable(model, vᵁ⁻[i=1:n] ≥ 0)
        # vᵢᴸ⁺ ≥ 0, vᵢᵁ⁺ ≥ 0
        @variable(model, vᴸ⁺[i=1:n] ≥ 0); @variable(model, vᵁ⁺[i=1:n] ≥ 0)

        for i = 1:n
            wᵢᴸ⁻ = wᴸ⁻[i]; wᵢᵁ⁻ = wᵁ⁻[i]
            wᵢᴸ⁺ = wᴸ⁺[i]; wᵢᵁ⁺ = wᵁ⁺[i]
            εᵢᴸ = εᴸ[i]; εᵢᵁ = εᵁ[i]
            vᵢᴸ⁻ = vᴸ⁻[i]; vᵢᵁ⁻ = vᵁ⁻[i]
            vᵢᴸ⁺ = vᴸ⁺[i]; vᵢᵁ⁺ = vᵁ⁺[i]

            @constraint(model, wᵢᵁ⁺ ≥ wᵢᵁ⁻)
            @constraint(model, wᵢᵁ⁻ ≥ wᵢᴸ⁻)
            @constraint(model, wᵢᴸ⁻ ≥ wᵢᴸ⁺)

            # 正規性条件
            ∑wⱼᴸ⁻ = sum(map(j -> wᴸ⁻[j] + εᴸ[j], filter(j -> i != j, 1:n)))
            @constraint(model, ∑wⱼᴸ⁻ + wᵢᵁ⁻ - εᵢᵁ ≤ 1)
            ∑wⱼᵁ⁻ = sum(map(j -> wᵁ⁻[j] - εᵁ[j], filter(j -> i != j, 1:n)))
            @constraint(model, ∑wⱼᵁ⁻ + wᵢᴸ⁻ + εᵢᴸ ≥ 1)
            ∑wⱼᴸ⁺ = sum(map(j -> wᴸ⁺[j] - εᴸ[j], filter(j -> i != j, 1:n)))
            @constraint(model, ∑wⱼᴸ⁺ + wᵢᵁ⁺ + εᵢᵁ ≤ 1)
            ∑wⱼᵁ⁺ = sum(map(j -> wᵁ⁺[j] + εᵁ[j], filter(j -> i != j, 1:n)))
            @constraint(model, ∑wⱼᵁ⁺ + wᵢᴸ⁺ - εᵢᴸ ≥ 1)

            for j = 1:n, k = 1:m
                if i == j continue end

                Aₖ = pcms[k]
                aₖᵢⱼᴸ = Aₖ[i,j].lo; aₖᵢⱼᵁ = Aₖ[i,j].hi
                wⱼᴸ⁻ = wᴸ⁻[j]; wⱼᵁ⁻ = wᵁ⁻[j]
                wⱼᴸ⁺ = wᴸ⁺[j]; wⱼᵁ⁺ = wᵁ⁺[j]

                @constraint(model, wᵢᴸ⁻ ≥ aₖᵢⱼᴸ * wⱼᵁ⁻ - εᵢᴸ)
                @constraint(model, wᵢᵁ⁻ ≤ aₖᵢⱼᵁ * wⱼᴸ⁻ + εᵢᵁ)
                @constraint(model, wᵢᴸ⁺ ≤ aₖᵢⱼᴸ * wⱼᵁ⁺ + εᵢᴸ)
                @constraint(model, wᵢᵁ⁺ ≥ aₖᵢⱼᵁ * wⱼᴸ⁺ - εᵢᵁ)

                @constraint(model, vᵢᴸ⁻ ≤ aₖᵢⱼᴸ * wⱼᵁ⁻)
                @constraint(model, vᵢᵁ⁻ ≥ aₖᵢⱼᵁ * wⱼᴸ⁻)
                @constraint(model, vᵢᴸ⁺ ≥ aₖᵢⱼᴸ * wⱼᵁ⁺)
                @constraint(model, vᵢᵁ⁺ ≤ aₖᵢⱼᵁ * wⱼᴸ⁺)
            end
        end
        # 中心総和 = 1
        @constraint(model, sum(wᴸ⁻) + sum(wᵁ⁻) == 2)
        @constraint(model, sum(wᴸ⁺) + sum(wᵁ⁺) == 2)

        # 目的関数 ∑(εᵢᴸ + εᵢᵁ)
        @objective(model, Min, sum(εᴸ) + sum(εᵁ))

        optimize!(model)

        optimalValue = sum(value.(εᴸ)) + sum(value.(εᵁ))

        return (
            wᴸ⁻=map(i -> correctPrecisionLoss(value(wᴸ⁻[i]), value(wᵁ⁻[i])), 1:n),
            wᵁ⁻=value.(wᵁ⁻),
            wᴸ⁺=map(i -> correctPrecisionLoss(value(wᴸ⁺[i]), value(wᵁ⁺[i])), 1:n),
            wᵁ⁺=value.(wᵁ⁺),
            εᴸ=value.(εᴸ), εᵁ=value.(εᵁ),
            vᴸ⁻=map(i -> correctPrecisionLoss(value(vᴸ⁻[i]), value(vᵁ⁻[i])), 1:n),
            vᵁ⁻=value.(vᵁ⁻),
            vᴸ⁺=value.(vᴸ⁺), vᵁ⁺=value.(vᴸ⁺),
            optimalValue=optimalValue
        )
    finally
        # エラー終了時にも変数などを消去する
        empty!(model)
    end
end

function importance2TwofoldIntervalPCM_m2(
        lpResult::ApproximationLPResult_m2{T}
        )::Matrix{TwofoldInterval{T}} where {T <: Real}
    n = length(lpResult.wᴸ⁻)

    wᴸ⁻ = lpResult.wᴸ⁻; wᵁ⁻ = lpResult.wᵁ⁻
    wᴸ⁺ = lpResult.wᴸ⁺; wᵁ⁺ = lpResult.wᵁ⁺

    A = fill((1..1, 1..1), (n,n))
    for i = 1:n, j = 1:n
        if i == j continue end

        wᵢᴸ⁻ = wᴸ⁻[i]; wᵢᵁ⁻ = wᵁ⁻[i]
        wᵢᴸ⁺ = wᴸ⁺[i]; wᵢᵁ⁺ = wᵁ⁺[i]
        wⱼᴸ⁻ = wᴸ⁻[j]; wⱼᵁ⁻ = wᵁ⁻[j]
        wⱼᴸ⁺ = wᴸ⁺[j]; wⱼᵁ⁺ = wᵁ⁺[j]

        aᵢⱼᴸ⁻ = wᵢᴸ⁻/wⱼᵁ⁻; aᵢⱼᵁ⁻ = wᵢᵁ⁻/wⱼᴸ⁻
        aᵢⱼᴸ⁺ = wᵢᴸ⁺/wⱼᵁ⁺; aᵢⱼᵁ⁺ = wᵢᵁ⁺/wⱼᴸ⁺

        aᵢⱼᴸ⁻ = correctPrecisionLoss(aᵢⱼᴸ⁻, aᵢⱼᴸ⁺)
        aᵢⱼᵁ⁻ = correctPrecisionLoss(aᵢⱼᵁ⁻, aᵢⱼᴸ⁻)
        aᵢⱼᵁ⁺ = correctPrecisionLoss(aᵢⱼᵁ⁺, aᵢⱼᵁ⁻)

        A[i,j] = (aᵢⱼᴸ⁻..aᵢⱼᵁ⁻, aᵢⱼᴸ⁺..aᵢⱼᵁ⁺)
    end

    if !isTwofoldIntervalPCM(A)
        throw(ArgumentError("Given matrix is not valid as twofold interval matrix."))
    end

    return A
end
