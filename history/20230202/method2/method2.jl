using JuMP
import HiGHS

include("../method1/index.jl")
include("../nearlyEqual/index.jl")
include("../twofoldInterval/index.jl")
include("../twofoldIntervalPCM/index.jl")

CancellingLP_m2 = LPResult_m1
solveCancellingLP_m2 = solveLP_m1
updatePCM_m2 = updatePCM_m1

function solveIndividualCancellingLP_m2(
        pcms::Vector{Matrix{Interval{T}}}
        )::T where {T <: Real}
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
        @variable(model, εᴸ⁻[k=1:m,i=1:n] ≥ 0) # εᵢᴸ⁻ ≥ 0
        @variable(model, εᵁ⁻[k=1:m,i=1:n] ≥ 0) # εᵢᵁ⁻ ≥ 0
        @variable(model, εᴸ⁺[k=1:m,i=1:n] ≥ 0) # εᵢᴸ⁺ ≥ 0
        @variable(model, εᵁ⁺[k=1:m,i=1:n] ≥ 0) # εᵢᵁ⁺ ≥ 0

        for i = 1:n
            wᵢᴸ⁻ = wᴸ⁻[i]; wᵢᵁ⁻ = wᵁ⁻[i]
            wᵢᴸ⁺ = wᴸ⁺[i]; wᵢᵁ⁺ = wᵁ⁺[i]

            @constraint(model, wᵢᵁ⁺ ≥ wᵢᵁ⁻)
            @constraint(model, wᵢᵁ⁻ ≥ wᵢᴸ⁻)
            @constraint(model, wᵢᴸ⁻ ≥ wᵢᴸ⁺)

            # 正規性条件
            ∑wⱼᴸ⁻ = sum(map(j -> wᴸ⁻[j], filter(j -> i != j, 1:n)))
            @constraint(model, ∑wⱼᴸ⁻ + wᵢᵁ⁻ ≤ 1)
            ∑wⱼᵁ⁻ = sum(map(j -> wᵁ⁻[j], filter(j -> i != j, 1:n)))
            @constraint(model, ∑wⱼᵁ⁻ + wᵢᴸ⁻ ≥ 1)
            ∑wⱼᴸ⁺ = sum(map(j -> wᴸ⁺[j], filter(j -> i != j, 1:n)))
            @constraint(model, ∑wⱼᴸ⁺ + wᵢᵁ⁺ ≤ 1)
            ∑wⱼᵁ⁺ = sum(map(j -> wᵁ⁺[j], filter(j -> i != j, 1:n)))
            @constraint(model, ∑wⱼᵁ⁺ + wᵢᴸ⁺ ≥ 1)

            for j = 1:n, k = 1:m
                if i == j continue end

                Aₖ = pcms[k]
                aₖᵢⱼᴸ = Aₖ[i,j].lo; aₖᵢⱼᵁ = Aₖ[i,j].hi
                wⱼᴸ⁻ = wᴸ⁻[j]; wⱼᵁ⁻ = wᵁ⁻[j]
                wⱼᴸ⁺ = wᴸ⁺[j]; wⱼᵁ⁺ = wᵁ⁺[j]
                εₖᵢᴸ⁻ = εᴸ⁻[k,i]; εₖᵢᵁ⁻ = εᵁ⁻[k,i]
                εₖᵢᴸ⁺ = εᴸ⁺[k,i]; εₖᵢᵁ⁺ = εᵁ⁺[k,i]

                @constraint(model, wᵢᴸ⁻ ≥ aₖᵢⱼᴸ * wⱼᵁ⁻ - εₖᵢᴸ⁻)
                @constraint(model, wᵢᵁ⁻ ≤ aₖᵢⱼᵁ * wⱼᴸ⁻ + εₖᵢᵁ⁻)
                @constraint(model, wᵢᴸ⁺ ≤ aₖᵢⱼᴸ * wⱼᵁ⁺ + εₖᵢᴸ⁺)
                @constraint(model, wᵢᵁ⁺ ≥ aₖᵢⱼᵁ * wⱼᴸ⁺ - εₖᵢᵁ⁺)
            end
        end
        # 中心総和 = 1
        @constraint(model, sum(wᴸ⁻) + sum(wᵁ⁻) == 2)
        @constraint(model, sum(wᴸ⁺) + sum(wᵁ⁺) == 2)

        # 目的関数 ∑(εₖᵢᴸ⁻ + εₖᵢᵁ⁻ + εₖᵢᴸ⁺ + εₖᵢᵁ⁺)
        @objective(model, Min, sum(εᴸ⁻) + sum(εᵁ⁻) + sum(εᴸ⁺) + sum(εᵁ⁺))

        optimize!(model)

        optimalValue =
            sum(value.(εᴸ⁻)) + sum(value.(εᵁ⁻)) +
            sum(value.(εᴸ⁺)) + sum(value.(εᵁ⁺))

        return optimalValue
    finally
        # エラー終了時にも変数などを消去する
        empty!(model)
    end
end

ApproximationLPResult_m2 = @NamedTuple{
    wᴸ⁻::Vector{T}, wᵁ⁻::Vector{T},
    wᴸ⁺::Vector{T}, wᵁ⁺::Vector{T},
    εᴸ⁻::Matrix{T}, εᵁ⁻::Matrix{T},
    εᴸ⁺::Matrix{T}, εᵁ⁺::Matrix{T},
    optimalValue::T
    } where {T <: Real}

function solveApproximationLP_m2(
        pcms::Vector{Matrix{Interval{T}}},
        εᶿ::T # 上付き * がないので θ で代用
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
        @variable(model, εᴸ⁻[k=1:m,i=1:n] ≥ 0) # εᵢᴸ⁻ ≥ 0
        @variable(model, εᵁ⁻[k=1:m,i=1:n] ≥ 0) # εᵢᵁ⁻ ≥ 0
        @variable(model, εᴸ⁺[k=1:m,i=1:n] ≥ 0) # εᵢᴸ⁺ ≥ 0
        @variable(model, εᵁ⁺[k=1:m,i=1:n] ≥ 0) # εᵢᵁ⁺ ≥ 0

        for i = 1:n
            wᵢᴸ⁻ = wᴸ⁻[i]; wᵢᵁ⁻ = wᵁ⁻[i]
            wᵢᴸ⁺ = wᴸ⁺[i]; wᵢᵁ⁺ = wᵁ⁺[i]

            @constraint(model, wᵢᵁ⁺ ≥ wᵢᵁ⁻)
            @constraint(model, wᵢᵁ⁻ ≥ wᵢᴸ⁻)
            @constraint(model, wᵢᴸ⁻ ≥ wᵢᴸ⁺)

            # 正規性条件
            ∑wⱼᴸ⁻ = sum(map(j -> wᴸ⁻[j], filter(j -> i != j, 1:n)))
            @constraint(model, ∑wⱼᴸ⁻ + wᵢᵁ⁻ ≤ 1)
            ∑wⱼᵁ⁻ = sum(map(j -> wᵁ⁻[j], filter(j -> i != j, 1:n)))
            @constraint(model, ∑wⱼᵁ⁻ + wᵢᴸ⁻ ≥ 1)
            ∑wⱼᴸ⁺ = sum(map(j -> wᴸ⁺[j], filter(j -> i != j, 1:n)))
            @constraint(model, ∑wⱼᴸ⁺ + wᵢᵁ⁺ ≤ 1)
            ∑wⱼᵁ⁺ = sum(map(j -> wᵁ⁺[j], filter(j -> i != j, 1:n)))
            @constraint(model, ∑wⱼᵁ⁺ + wᵢᴸ⁺ ≥ 1)

            for j = 1:n, k = 1:m
                if i == j continue end

                Aₖ = pcms[k]
                aₖᵢⱼᴸ = Aₖ[i,j].lo; aₖᵢⱼᵁ = Aₖ[i,j].hi
                wⱼᴸ⁻ = wᴸ⁻[j]; wⱼᵁ⁻ = wᵁ⁻[j]
                wⱼᴸ⁺ = wᴸ⁺[j]; wⱼᵁ⁺ = wᵁ⁺[j]
                εₖᵢᴸ⁻ = εᴸ⁻[k,i]; εₖᵢᵁ⁻ = εᵁ⁻[k,i]
                εₖᵢᴸ⁺ = εᴸ⁺[k,i]; εₖᵢᵁ⁺ = εᵁ⁺[k,i]

                @constraint(model, wᵢᴸ⁻ ≥ aₖᵢⱼᴸ * wⱼᵁ⁻ - εₖᵢᴸ⁻)
                @constraint(model, wᵢᵁ⁻ ≤ aₖᵢⱼᵁ * wⱼᴸ⁻ + εₖᵢᵁ⁻)
                @constraint(model, wᵢᴸ⁺ ≤ aₖᵢⱼᴸ * wⱼᵁ⁺ + εₖᵢᴸ⁺)
                @constraint(model, wᵢᵁ⁺ ≥ aₖᵢⱼᵁ * wⱼᴸ⁺ - εₖᵢᵁ⁺)
            end
        end
        # 中心総和 = 1
        @constraint(model, sum(wᴸ⁻) + sum(wᵁ⁻) == 2)
        @constraint(model, sum(wᴸ⁺) + sum(wᵁ⁺) == 2)
        # 修正量の総和の制約
        @constraint(model, sum(εᴸ⁻) + sum(εᵁ⁻) + sum(εᴸ⁺) + sum(εᵁ⁺) == εᶿ)

        # 目的関数 ∑(wᵢᴸ⁻ - wᵢᴸ⁺ + wᵢᵁ⁺ - wᵢᵁ⁻)
        @objective(model, Min, sum(wᴸ⁻) - sum(wᴸ⁺) + sum(wᵁ⁺) - sum(wᵁ⁻))

        optimize!(model)

        optimalValue =
            sum(value.(wᴸ⁻)) - sum(value.(wᴸ⁺)) +
            sum(value.(wᵁ⁺)) - sum(value.(wᵁ⁻))

        return (
            wᴸ⁻=map(i -> correctPrecisionLoss(value(wᴸ⁻[i]), value(wᵁ⁻[i])), 1:n),
            wᵁ⁻=value.(wᵁ⁻),
            wᴸ⁺=map(i -> correctPrecisionLoss(value(wᴸ⁺[i]), value(wᵁ⁺[i])), 1:n),
            wᵁ⁺=value.(wᵁ⁺),
            εᴸ⁻=value.(εᴸ⁻), εᵁ⁻=value.(εᵁ⁻),
            εᴸ⁺=value.(εᴸ⁺), εᵁ⁺=value.(εᵁ⁺),
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

        if aᵢⱼᴸ⁻ > aᵢⱼᵁ⁻
            A[i, j] = (emptyinterval(), aᵢⱼᴸ⁺..aᵢⱼᵁ⁺)
        else
            A[i, j] = (aᵢⱼᴸ⁻..aᵢⱼᵁ⁻, aᵢⱼᴸ⁺..aᵢⱼᵁ⁺)
        end
    end

    if !isTwofoldIntervalPCM(A)
        throw(ArgumentError("Given matrix is not valid as twofold interval matrix."))
    end

    return A
end
