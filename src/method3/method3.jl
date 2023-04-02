using IntervalArithmetic
using JuMP
import HiGHS

include("../nearlyEqual/index.jl")
include("../twofoldInterval/index.jl")
include("../twofoldIntervalPCM/index.jl")

function solveFeasibilityCheckLP_m3(
        Aₖ::Matrix{Interval{T}})::T where {T <: Real}
    ε = 1e-8

    if !isIntervalPCM(Aₖ)
        throw(ArgumentError("Given matrix is not valid as interval matrix."))
    end

    m, n = size(Aₖ)
    model = Model(HiGHS.Optimizer)
    set_silent(model)

    try
        # wₖᵢᴸ⁻ ≥ ε, wₖᵢᵁ⁻ ≥ ε
        @variable(model, wₖᴸ⁻[i=1:n] ≥ ε); @variable(model, wₖᵁ⁻[i=1:n] ≥ ε)
        # wₖᵢᴸ⁺ ≥ ε, wₖᵢᵁ⁺ ≥ ε
        @variable(model, wₖᴸ⁺[i=1:n] ≥ ε); @variable(model, wₖᵁ⁺[i=1:n] ≥ ε)
        # δₖᵢᴸ ≥ 0, δₖᵢᵁ ≥ 0
        @variable(model, δₖᴸ[i=1:n] ≥ 0); @variable(model, δₖᵁ[i=1:n] ≥ 0)

        for i = 1:n
            wₖᵢᴸ⁻ = wₖᴸ⁻[i]; wₖᵢᵁ⁻ = wₖᵁ⁻[i]; wₖᵢᴸ⁺ = wₖᴸ⁺[i]; wₖᵢᵁ⁺ = wₖᵁ⁺[i]
            δₖᵢᴸ = δₖᴸ[i]; δₖᵢᵁ = δₖᵁ[i]

            @constraint(model, wₖᵢᵁ⁺ ≥ wₖᵢᵁ⁻)
            @constraint(model, wₖᵢᵁ⁻ ≥ wₖᵢᴸ⁻)
            @constraint(model, wₖᵢᴸ⁻ ≥ wₖᵢᴸ⁺)

            # 正規性条件
            ∑wₖⱼᴸ⁻ = sum(map(j -> wₖᴸ⁻[j], filter(j -> i != j, 1:n)))
            @constraint(model, ∑wₖⱼᴸ⁻ + wₖᵢᵁ⁻ ≤ 1)
            ∑wₖⱼᵁ⁻ = sum(map(j -> wₖᵁ⁻[j], filter(j -> i != j, 1:n)))
            @constraint(model, ∑wₖⱼᵁ⁻ + wₖᵢᴸ⁻ ≥ 1)
            ∑wₖⱼᴸ⁺ = sum(map(j -> wₖᴸ⁺[j], filter(j -> i != j, 1:n)))
            @constraint(model, ∑wₖⱼᴸ⁺ + wₖᵢᵁ⁺ ≤ 1)
            ∑wₖⱼᵁ⁺ = sum(map(j -> wₖᵁ⁺[j], filter(j -> i != j, 1:n)))
            @constraint(model, ∑wₖⱼᵁ⁺ + wₖᵢᴸ⁺ ≥ 1)

            for j = 1:n
                if i == j continue end

                aₖᵢⱼᴸ = Aₖ[i,j].lo; aₖᵢⱼᵁ = Aₖ[i,j].hi
                wₖⱼᴸ⁻ = wₖᴸ⁻[j]; wₖⱼᵁ⁻ = wₖᵁ⁻[j]; wₖⱼᴸ⁺ = wₖᴸ⁺[j]; wₖⱼᵁ⁺ = wₖᵁ⁺[j]

                @constraint(model, wₖᵢᴸ⁺ ≤ aₖᵢⱼᴸ * wₖⱼᵁ⁺)
                @constraint(model, wₖᵢᴸ⁻ ≥ aₖᵢⱼᴸ * wₖⱼᵁ⁻ - δₖᵢᴸ)
                @constraint(model, wₖᵢᵁ⁻ ≤ aₖᵢⱼᵁ * wₖⱼᴸ⁻ + δₖᵢᵁ)
                @constraint(model, wₖᵢᵁ⁺ ≥ aₖᵢⱼᵁ * wₖⱼᴸ⁺)
            end
        end
        # 中心総和 = 1
        @constraint(model, sum(wₖᴸ⁻) + sum(wₖᵁ⁻) == 2)
        @constraint(model, sum(wₖᴸ⁺) + sum(wₖᵁ⁺) == 2)

        # 目的関数 ∑(δₖᵢᴸ + δₖᵢᵁ)
        @objective(model, Min, sum(δₖᴸ) + sum(δₖᵁ))

        optimize!(model)

        optimalValue = sum(value.(δₖᴸ)) + sum(value.(δₖᵁ))

        return optimalValue
    finally
        # エラー終了時にも変数などを消去する
        empty!(model)
    end
end

ApproximationLPResult_m3 = @NamedTuple{
    wₖᴸ⁻::Vector{T}, wₖᵁ⁻::Vector{T},
    wₖᴸ⁺::Vector{T}, wₖᵁ⁺::Vector{T},
    optimalValue::T
    } where {T <: Real}

function solveUpperApproximationLP_m3(
        Aₖ::Matrix{Interval{T}}
        )::ApproximationLPResult_m3{T} where {T <: Real}
    ε = 1e-8

    if !isIntervalPCM(Aₖ)
        throw(ArgumentError("Given matrix is not valid as interval matrix."))
    end

    m, n = size(Aₖ)
    model = Model(HiGHS.Optimizer)
    set_silent(model)

    try
        # wₖᵢᴸ⁺ ≥ ε, wₖᵢᵁ⁺ ≥ ε
        @variable(model, wₖᴸ⁺[i=1:n] ≥ ε); @variable(model, wₖᵁ⁺[i=1:n] ≥ ε)

        for i = 1:n
            wₖᵢᴸ⁺ = wₖᴸ⁺[i]; wₖᵢᵁ⁺ = wₖᵁ⁺[i]

            @constraint(model, wₖᵢᵁ⁺ ≥ wₖᵢᴸ⁺)

            # 正規性条件
            ∑wₖⱼᴸ⁺ = sum(map(j -> wₖᴸ⁺[j], filter(j -> i != j, 1:n)))
            @constraint(model, ∑wₖⱼᴸ⁺ + wₖᵢᵁ⁺ ≤ 1)
            ∑wₖⱼᵁ⁺ = sum(map(j -> wₖᵁ⁺[j], filter(j -> i != j, 1:n)))
            @constraint(model, ∑wₖⱼᵁ⁺ + wₖᵢᴸ⁺ ≥ 1)

            for j = 1:n
                if i == j continue end

                aₖᵢⱼᴸ = Aₖ[i,j].lo; aₖᵢⱼᵁ = Aₖ[i,j].hi
                wₖⱼᴸ⁺ = wₖᴸ⁺[j]; wₖⱼᵁ⁺ = wₖᵁ⁺[j]

                @constraint(model, wₖᵢᴸ⁺ ≤ aₖᵢⱼᴸ * wₖⱼᵁ⁺)
                @constraint(model, wₖᵢᵁ⁺ ≥ aₖᵢⱼᵁ * wₖⱼᴸ⁺)
            end
        end
        @constraint(model, sum(wₖᴸ⁺) + sum(wₖᵁ⁺) == 2) # 中心総和 = 1

        # 目的関数 ∑(wₖᵢᴸ⁻ - wₖᵢᴸ⁺ + wₖᵢᵁ⁺ - wₖᵢᵁ⁻)
        @objective(model, Min, -sum(wₖᴸ⁺) + sum(wₖᵁ⁺))

        optimize!(model)

        optimalValue = -sum(value.(wₖᴸ⁺)) + sum(value.(wₖᵁ⁺))

        return (
            # IntervalArithmetic の ∅ に合わせて
            # wₖᵢᴸ⁻ = ∞, wₖᵢᵁ⁻ = -∞
            wₖᴸ⁻=fill(∞, n), wₖᵁ⁻=fill(-∞, n),
            wₖᴸ⁺=map(i -> correctPrecisionLoss(value(wₖᴸ⁺[i]), value(wₖᵁ⁺[i])), 1:n),
            wₖᵁ⁺=value.(wₖᵁ⁺),
            optimalValue=optimalValue
        )
    finally
        # エラー終了時にも変数などを消去する
        empty!(model)
    end
end

function solveBothApproximationLP_m3(
        Aₖ::Matrix{Interval{T}}
        )::ApproximationLPResult_m3{T} where {T <: Real}
    ε = 1e-8

    if !isIntervalPCM(Aₖ)
        throw(ArgumentError("Given matrix is not valid as interval matrix."))
    end

    m, n = size(Aₖ)
    model = Model(HiGHS.Optimizer)
    set_silent(model)

    try
        # wₖᵢᴸ⁻ ≥ ε, wₖᵢᵁ⁻ ≥ ε
        @variable(model, wₖᴸ⁻[i=1:n] ≥ ε); @variable(model, wₖᵁ⁻[i=1:n] ≥ ε)
        # wₖᵢᴸ⁺ ≥ ε, wₖᵢᵁ⁺ ≥ ε
        @variable(model, wₖᴸ⁺[i=1:n] ≥ ε); @variable(model, wₖᵁ⁺[i=1:n] ≥ ε)

        for i = 1:n
            wₖᵢᴸ⁻ = wₖᴸ⁻[i]; wₖᵢᵁ⁻ = wₖᵁ⁻[i]; wₖᵢᴸ⁺ = wₖᴸ⁺[i]; wₖᵢᵁ⁺ = wₖᵁ⁺[i]

            @constraint(model, wₖᵢᵁ⁺ ≥ wₖᵢᵁ⁻)
            @constraint(model, wₖᵢᵁ⁻ ≥ wₖᵢᴸ⁻)
            @constraint(model, wₖᵢᴸ⁻ ≥ wₖᵢᴸ⁺)

            # 正規性条件
            ∑wₖⱼᴸ⁻ = sum(map(j -> wₖᴸ⁻[j], filter(j -> i != j, 1:n)))
            @constraint(model, ∑wₖⱼᴸ⁻ + wₖᵢᵁ⁻ ≤ 1)
            ∑wₖⱼᵁ⁻ = sum(map(j -> wₖᵁ⁻[j], filter(j -> i != j, 1:n)))
            @constraint(model, ∑wₖⱼᵁ⁻ + wₖᵢᴸ⁻ ≥ 1)
            ∑wₖⱼᴸ⁺ = sum(map(j -> wₖᴸ⁺[j], filter(j -> i != j, 1:n)))
            @constraint(model, ∑wₖⱼᴸ⁺ + wₖᵢᵁ⁺ ≤ 1)
            ∑wₖⱼᵁ⁺ = sum(map(j -> wₖᵁ⁺[j], filter(j -> i != j, 1:n)))
            @constraint(model, ∑wₖⱼᵁ⁺ + wₖᵢᴸ⁺ ≥ 1)

            for j = 1:n
                if i == j continue end

                aₖᵢⱼᴸ = Aₖ[i,j].lo; aₖᵢⱼᵁ = Aₖ[i,j].hi
                wₖⱼᴸ⁻ = wₖᴸ⁻[j]; wₖⱼᵁ⁻ = wₖᵁ⁻[j]; wₖⱼᴸ⁺ = wₖᴸ⁺[j]; wₖⱼᵁ⁺ = wₖᵁ⁺[j]

                @constraint(model, wₖᵢᴸ⁺ ≤ aₖᵢⱼᴸ * wₖⱼᵁ⁺)
                @constraint(model, wₖᵢᴸ⁻ ≥ aₖᵢⱼᴸ * wₖⱼᵁ⁻)
                @constraint(model, wₖᵢᵁ⁻ ≤ aₖᵢⱼᵁ * wₖⱼᴸ⁻)
                @constraint(model, wₖᵢᵁ⁺ ≥ aₖᵢⱼᵁ * wₖⱼᴸ⁺)
            end
        end
        # 中心総和 = 1
        @constraint(model, sum(wₖᴸ⁻) + sum(wₖᵁ⁻) == 2)
        @constraint(model, sum(wₖᴸ⁺) + sum(wₖᵁ⁺) == 2)

        # 目的関数 ∑(wₖᵢᴸ⁻ - wₖᵢᴸ⁺ + wₖᵢᵁ⁺ - wₖᵢᵁ⁻)
        @objective(model, Min, sum(wₖᴸ⁻) - sum(wₖᴸ⁺) + sum(wₖᵁ⁺) - sum(wₖᵁ⁻))

        optimize!(model)

        optimalValue = sum(value.(wₖᴸ⁻)) - sum(value.(wₖᴸ⁺)) + sum(value.(wₖᵁ⁺)) - sum(value.(wₖᵁ⁻))

        return (
            wₖᴸ⁻=map(i -> correctPrecisionLoss(value(wₖᴸ⁻[i]), value(wₖᵁ⁻[i])), 1:n),
            wₖᵁ⁻=value.(wₖᵁ⁻),
            wₖᴸ⁺=map(i -> correctPrecisionLoss(value(wₖᴸ⁺[i]), value(wₖᵁ⁺[i])), 1:n),
            wₖᵁ⁺=value.(wₖᵁ⁺),
            optimalValue=optimalValue
        )
    finally
        # エラー終了時にも変数などを消去する
        empty!(model)
    end
end

function solveApproximationLP_m3(
        Aₖ::Matrix{Interval{T}}
        )::ApproximationLPResult_m3{T} where {T <: Real}
    # W⁻ = ∅ かどうか判定するLPの最適値が0ならば W⁻ ≠ ∅
    if solveFeasibilityCheckLP_m3(Aₖ) == 0
        return solveBothApproximationLP_m3(Aₖ)
    else
        return solveUpperApproximationLP_m3(Aₖ)
    end
end

TBoundaries_m3 = @NamedTuple{
    tₖᴸ⁻::T, tₖᵁ⁻::T, tₖᴸ⁺::T, tₖᵁ⁺::T,
    } where {T <: Real}

function calculateTBoundaries_m3(
        lpResult::ApproximationLPResult_m3{T}
        )::TBoundaries_m3{T} where {T <: Real}
    wₖᴸ⁻ = lpResult.wₖᴸ⁻; wₖᵁ⁻ = lpResult.wₖᵁ⁻
    wₖᴸ⁺ = lpResult.wₖᴸ⁺; wₖᵁ⁺ = lpResult.wₖᵁ⁺

    n = length(wₖᴸ⁻)

    # W⁻ = ∅
    if any(isinf.(wₖᴸ⁻)) || any(isinf.(wₖᵁ⁻))
        tₖᴸ⁺ = 1 / minimum(i -> sum(map(j -> wₖᵁ⁺[j], filter(j -> i != j, 1:n))) + wₖᴸ⁺[i], 1:n)
        tₖᵁ⁺ = 1 / maximum(i -> sum(map(j -> wₖᴸ⁺[j], filter(j -> i != j, 1:n))) + wₖᵁ⁺[i], 1:n)
    
        return (tₖᴸ⁻=-∞, tₖᵁ⁻=∞, tₖᴸ⁺=tₖᴸ⁺, tₖᵁ⁺=tₖᵁ⁺)
    end

    tₖᴸ⁻ = 1 / minimum(i -> sum(map(j -> wₖᵁ⁻[j], filter(j -> i != j, 1:n))) + wₖᴸ⁻[i], 1:n)
    tₖᵁ⁻ = 1 / maximum(i -> sum(map(j -> wₖᴸ⁻[j], filter(j -> i != j, 1:n))) + wₖᵁ⁻[i], 1:n)
    tₖᴸ⁺ = 1 / minimum(i -> sum(map(j -> wₖᵁ⁺[j], filter(j -> i != j, 1:n))) + wₖᴸ⁺[i], 1:n)
    tₖᵁ⁺ = 1 / maximum(i -> sum(map(j -> wₖᴸ⁺[j], filter(j -> i != j, 1:n))) + wₖᵁ⁺[i], 1:n)

    return (tₖᴸ⁻=tₖᴸ⁻, tₖᵁ⁻=tₖᵁ⁻, tₖᴸ⁺=tₖᴸ⁺, tₖᵁ⁺=tₖᵁ⁺)
end

ConcatLPResult_m3 = @NamedTuple{
    t⁻::Vector{T}, t⁺::Vector{T},
    wᴸ::Vector{T}, wᵁ::Vector{T},
    vᴸ⁻::Vector{T}, vᵁ⁻::Vector{T},
    vᴸ⁺::Vector{T}, vᵁ⁺::Vector{T},
    εᴸ::Vector{T}, εᵁ::Vector{T},
    optimalValue::T
    } where {T <: Real}

function solveConcatLP_m3(
        lpResults::AbstractArray{ApproximationLPResult_m3{T}, 1},
        tBoundaries::AbstractArray{TBoundaries_m3{T}, 1}
        )::ConcatLPResult_m3{T} where {T <: Real}
    ε = 1e-8

    m = length(lpResults)
    n = length(lpResults[1].wₖᴸ⁻)

    model = Model(HiGHS.Optimizer)
    set_silent(model)

    try
        # tₖ⁻ ≥ ε, tₖ⁺ ≥ ε
        @variable(model, t⁻[i=1:m] ≥ ε); @variable(model, t⁺[i=1:m] ≥ ε)
        # wᵢᴸ ≥ ε, wᵢᵁ ≥ ε
        @variable(model, wᴸ[i=1:n] ≥ ε); @variable(model, wᵁ[i=1:n] ≥ ε)
        # vᵢᴸ⁻ ≥ ε, vᵢᵁ⁻ ≥ ε
        @variable(model, vᴸ⁻[i=1:n] ≥ ε); @variable(model, vᵁ⁻[i=1:n] ≥ ε)
        # vᵢᴸ⁺ ≥ ε, vᵢᵁ⁺ ≥ ε
        @variable(model, vᴸ⁺[i=1:n] ≥ ε); @variable(model, vᵁ⁺[i=1:n] ≥ ε)
        # εᵢᴸ ≥ 0, εᵢᵁ ≥ 0
        @variable(model, εᴸ[i=1:n] ≥ 0); @variable(model, εᵁ[i=1:n] ≥ 0)

        for i = 1:n
            wᵢᴸ = wᴸ[i]; wᵢᵁ = wᵁ[i]
            vᵢᴸ⁻ = vᴸ⁻[i]; vᵢᵁ⁻ = vᵁ⁻[i]; vᵢᴸ⁺ = vᴸ⁺[i]; vᵢᵁ⁺ = vᵁ⁺[i]
            εᵢᴸ = εᴸ[i]; εᵢᵁ = εᵁ[i]

            # @constraint(model, vᵢᴸ⁻ ≤ maximum(filter(v -> v != 0 && isfinite(v), [
            #     map(k -> tBoundaries[k].tₖᵁ⁻ * lpResults[k].wₖᴸ⁻[i], 1:m)...,
            #     map(k -> tBoundaries[k].tₖᵁ⁺ * lpResults[k].wₖᴸ⁺[i], 1:m)...
            # ])))
            # @constraint(model, vᵢᵁ⁻ ≥ minimum(filter(v -> v != 0 && isfinite(v), [
            #     map(k -> tBoundaries[k].tₖᴸ⁻ * lpResults[k].wₖᵁ⁻[i], 1:m)...,
            #     map(k -> tBoundaries[k].tₖᴸ⁺ * lpResults[k].wₖᵁ⁺[i], 1:m)...
            # ])))

            @constraint(model, wᵢᵁ ≥ wᵢᴸ)
            @constraint(model, εᵢᴸ ≥ vᵢᴸ⁺ - vᵢᴸ⁻)
            @constraint(model, εᵢᵁ ≥ vᵢᵁ⁻ - vᵢᵁ⁺)

            # 正規性条件
            ∑wⱼᴸ = sum(map(j -> wᴸ[j], filter(j -> i != j, 1:n)))
            @constraint(model, ∑wⱼᴸ + wᵢᵁ ≤ 1)
            ∑wⱼᵁ = sum(map(j -> wᵁ[j], filter(j -> i != j, 1:n)))
            @constraint(model, ∑wⱼᵁ + wᵢᴸ ≥ 1)
        end
        @constraint(model, sum(wᴸ) + sum(wᵁ) == 2) # 中心総和 = 1

        for k = 1:m
            tₖ⁻ = t⁻[k]; tₖ⁺ = t⁺[k]
            tₖᴸ⁻ = tBoundaries[k].tₖᴸ⁻; tₖᵁ⁻ = tBoundaries[k].tₖᵁ⁻
            tₖᴸ⁺ = tBoundaries[k].tₖᴸ⁺; tₖᵁ⁺ = tBoundaries[k].tₖᵁ⁺

            # k ∈ M'
            if isfinite(tₖᴸ⁻) && isfinite(tₖᵁ⁻)
                @constraint(model, tₖᴸ⁻ ≤ tₖ⁻); @constraint(model, tₖ⁻ ≤ tₖᵁ⁻)
            end
            @constraint(model, tₖᴸ⁺ ≤ tₖ⁺); @constraint(model, tₖ⁺ ≤ tₖᵁ⁺)

            for i = 1:n
                wₖᵢᴸ⁻ = lpResults[k].wₖᴸ⁻[i]; wₖᵢᵁ⁻ = lpResults[k].wₖᵁ⁻[i]
                wₖᵢᴸ⁺ = lpResults[k].wₖᴸ⁺[i]; wₖᵢᵁ⁺ = lpResults[k].wₖᵁ⁺[i]
                vᵢᴸ⁻ = vᴸ⁻[i]; vᵢᵁ⁻ = vᵁ⁻[i]; vᵢᴸ⁺ = vᴸ⁺[i]; vᵢᵁ⁺ = vᵁ⁺[i]
                εᵢᴸ = εᴸ[i]; εᵢᵁ = εᵁ[i]
                wᵢᴸ = wᴸ[i]; wᵢᵁ = wᵁ[i]

                @constraint(model, tₖ⁺ * wₖᵢᴸ⁺ - εᵢᴸ ≤ wᵢᴸ)
                @constraint(model, wᵢᵁ ≤ tₖ⁺ * wₖᵢᵁ⁺ + εᵢᵁ)

                @constraint(model, vᵢᴸ⁺ ≥ tₖ⁺ * wₖᵢᴸ⁺)
                @constraint(model, vᵢᵁ⁺ ≤ tₖ⁺ * wₖᵢᵁ⁺)

                # k ∈ M'
                if isfinite(tₖᴸ⁻) && isfinite(tₖᵁ⁻)
                    @constraint(model, wᵢᴸ ≤ tₖ⁻ * wₖᵢᴸ⁻ + εᵢᴸ)
                    @constraint(model, tₖ⁻ * wₖᵢᵁ⁻ - εᵢᵁ ≤ wᵢᵁ)

                    @constraint(model, vᵢᴸ⁻ ≤ tₖ⁻ * wₖᵢᴸ⁻)
                    @constraint(model, vᵢᵁ⁻ ≥ tₖ⁻ * wₖᵢᵁ⁻)
                end
            end
        end

        # 目的関数 ∑(εᵢᴸ + εᵢᵁ)
        @objective(model, Min, sum(εᴸ) + sum(εᵁ))

        optimize!(model)

        optimalValue = sum(value.(εᴸ)) + sum(value.(εᵁ))

        return (
            t⁻=value.(t⁻), t⁺=value.(t⁺),
            # precision error の補正
            # wᵢᴸ と wᵢᵁ が十分に近い値ならば wᵢᴸ <- wᵢᵁ
            wᴸ=map(i -> correctPrecisionLoss(value(wᴸ[i]), value(wᵁ[i])), 1:n),
            wᵁ=value.(wᵁ),
            # precision error の補正
            # vᵢᴸ⁻ と vᵢᵁ⁻ が十t⁻分に近い値ならば vᵢᴸ⁻ <- vᵢᵁ⁻
            vᴸ⁻=map(i -> correctPrecisionLoss(value(vᴸ⁻[i]), value(vᵁ⁻[i])), 1:n),
            vᵁ⁻=value.(vᵁ⁻),
            vᴸ⁺=map(i -> correctPrecisionLoss(value(vᴸ⁺[i]), value(vᵁ⁺[i])), 1:n),
            vᵁ⁺=value.(vᵁ⁺),
            εᴸ=value.(εᴸ), εᵁ=value.(εᵁ),
            optimalValue=optimalValue
        )
    finally
        # エラー終了時にも変数などを消去する
        empty!(model)
    end
end

function generatePCM_m3(
        lpResult::ConcatLPResult_m3{T}
        )::Matrix{TwofoldInterval{T}} where {T <: Real}
    n = length(lpResult.wᴸ)

    wᴸ = lpResult.wᴸ; wᵁ = lpResult.wᵁ
    vᴸ⁻ = lpResult.vᴸ⁻; vᵁ⁻ = lpResult.vᵁ⁻
    vᴸ⁺ = lpResult.vᴸ⁺; vᵁ⁺ = lpResult.vᵁ⁺
    εᴸ = lpResult.εᴸ; εᵁ = lpResult.εᵁ

    # ([1, 1], [1, 1])で埋める
    Â = [(1..1, 1..1) for a in 1:n, b in 1:n]

    for i = 1:n, j = 1:n
        # 対角成分は (1..1, 1..1) で固定なので更新不要
        if i == j continue end

        wᵢᴸ = wᴸ[i]; wᵢᵁ = wᵁ[i]
        wⱼᴸ = wᴸ[j]; wⱼᵁ = wᵁ[j]

        vᵢᴸ⁻ = vᴸ⁻[i]; vᵢᵁ⁻ = vᵁ⁻[i]; vᵢᴸ⁺ = vᴸ⁺[i]; vᵢᵁ⁺ = vᵁ⁺[i]
        vⱼᴸ⁻ = vᴸ⁻[j]; vⱼᵁ⁻ = vᵁ⁻[j]; vⱼᴸ⁺ = vᴸ⁺[j]; vⱼᵁ⁺ = vᵁ⁺[j]
        εᵢᴸ = εᴸ[i]; εᵢᵁ = εᵁ[i]; εⱼᴸ = εᴸ[j]; εⱼᵁ = εᵁ[j]

        âᵢⱼᴸ⁻ = max((vᵢᴸ⁻ + εᵢᴸ) / wⱼᵁ, wᵢᴸ / (vⱼᵁ⁻ - εⱼᵁ))
        âᵢⱼᵁ⁻ = min((vᵢᵁ⁻ - εᵢᵁ) / wⱼᴸ, wᵢᵁ / (vⱼᴸ⁻ + εⱼᴸ))
        âᵢⱼᴸ⁺ = min((vᵢᴸ⁺ - εᵢᴸ) / wⱼᵁ, wᵢᴸ / (vⱼᵁ⁺ + εⱼᵁ))
        âᵢⱼᵁ⁺ = max((vᵢᵁ⁺ + εᵢᵁ) / wⱼᴸ, wᵢᵁ / (vⱼᴸ⁺ - εⱼᴸ))

        # âᵢⱼᵁ⁺ = âᵢⱼᴸ⁺ の場合などに precision error で âᵢⱼᵁ⁺ < âᵢⱼᴸ⁺ となることがある
        âᵢⱼᴸ⁻ = correctPrecisionLoss(âᵢⱼᴸ⁻, âᵢⱼᴸ⁺)
        âᵢⱼᵁ⁻ = correctPrecisionLoss(âᵢⱼᵁ⁻, âᵢⱼᴸ⁻)
        âᵢⱼᵁ⁺ = correctPrecisionLoss(âᵢⱼᵁ⁺, âᵢⱼᵁ⁻)
        
        âᵢⱼᵁ⁺ = correctPrecisionLoss(âᵢⱼᵁ⁺, âᵢⱼᴸ⁺)

        # (Âᵢⱼ⁻, Âᵢⱼ⁺)
        if âᵢⱼᴸ⁻ > âᵢⱼᵁ⁻
            Â[i, j] = (emptyinterval(), âᵢⱼᴸ⁺..âᵢⱼᵁ⁺)
        else
            Â[i, j] = (âᵢⱼᴸ⁻..âᵢⱼᵁ⁻, âᵢⱼᴸ⁺..âᵢⱼᵁ⁺)
        end
    end

    return Â
end
