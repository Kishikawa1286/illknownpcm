using IntervalArithmetic
using JuMP
import HiGHS

include("./solution.jl")
include("../twofoldInterval/index.jl")
include("../twofoldIntervalPCM/index.jl")

function solveLinear(
        A::Matrix{TwofoldInterval{T}},
        )::Solution{T} where {T <: Real}
    ε = 1e-8

    if !isTwofoldIntervalPCM(A)
        throw(ArgumentError("Given matrix is not valid as twofold interval matrix."))
    end

    m, n = size(A)
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    
    try
        # wᵢᴸ ≥ ε, wᵢᵁ ≥ ε
        @variable(model, wᴸ[i=1:n] ≥ ε); @variable(model, wᵁ[i=1:n] ≥ ε)
        # wᵢᴸ⁻ ≥ 0, wᵢᵁ⁻ ≥ 0
        @variable(model, wᴸ⁻[i=1:n] ≥ ε); @variable(model, wᵁ⁻[i=1:n] ≥ ε)
        # wᵢᴸ⁺ ≥ 0, wᵢᵁ⁺ ≥ 0
        @variable(model, wᴸ⁺[i=1:n] ≥ ε); @variable(model, wᵁ⁺[i=1:n] ≥ ε)
        # εᵢᴸ ≥ 0, εᵢᵁ ≥ 0
        @variable(model, εᴸ[i=1:n] ≥ 0); @variable(model, εᵁ[i=1:n] ≥ 0)

        for i = 1:n
            wᵢᴸ = wᴸ[i]; wᵢᵁ = wᵁ[i]
            wᵢᴸ⁻ = wᴸ⁻[i]; wᵢᵁ⁻ = wᵁ⁻[i]; wᵢᴸ⁺ = wᴸ⁺[i]; wᵢᵁ⁺ = wᵁ⁺[i]
            εᵢᴸ = εᴸ[i]; εᵢᵁ = εᵁ[i]

            @constraint(model, wᵢᵁ ≥ wᵢᴸ)
            @constraint(model, εᵢᴸ ≥ wᵢᴸ⁺ - wᵢᴸ⁻)
            @constraint(model, εᵢᵁ ≥ wᵢᵁ⁻ - wᵢᵁ⁺)

            # 正規性条件
            ∑wⱼᴸ = sum(map(j -> wᴸ[j], filter(j -> i != j, 1:n)))
            @constraint(model, ∑wⱼᴸ + wᵢᵁ ≤ 1)
            ∑wⱼᵁ = sum(map(j -> wᵁ[j], filter(j -> i != j, 1:n)))
            @constraint(model, ∑wⱼᵁ + wᵢᴸ ≥ 1)

            for j = 1:n
                if i == j
                    continue
                end

                aᵢⱼᴸ⁻ = A[i,j][1].lo; aᵢⱼᵁ⁻ = A[i,j][1].hi
                aᵢⱼᴸ⁺ = A[i,j][2].lo; aᵢⱼᵁ⁺ = A[i,j][2].hi
                wⱼᴸ = wᴸ[j]; wⱼᵁ = wᵁ[j]

                @constraint(model, aᵢⱼᴸ⁺ * wⱼᵁ - εᵢᴸ ≤ wᵢᴸ)
                @constraint(model, wᵢᴸ ≤ aᵢⱼᴸ⁻ * wⱼᵁ + εᵢᴸ)
                @constraint(model, aᵢⱼᵁ⁻ * wⱼᴸ - εᵢᵁ ≤ wᵢᵁ)
                @constraint(model, wᵢᵁ ≤ aᵢⱼᵁ⁺ * wⱼᴸ + εᵢᵁ)
                @constraint(model, wᵢᴸ⁻ ≤ aᵢⱼᴸ⁻ * wⱼᵁ)
                @constraint(model, wᵢᴸ⁺ ≥ aᵢⱼᴸ⁺ * wⱼᵁ)
                @constraint(model, wᵢᵁ⁻ ≥ aᵢⱼᵁ⁻ * wⱼᴸ)
                @constraint(model, wᵢᵁ⁺ ≤ aᵢⱼᵁ⁺ * wⱼᴸ)
            end
        end
        @constraint(model, sum(wᴸ) + sum(wᵁ) == 2) # 中心総和 = 1

        # 目的関数 ∑(εᵢᴸ + εᵢᵁ)
        @objective(model, Min, sum(εᴸ) + sum(εᵁ))

        JuMP.optimize!(model)

        optimalValue = sum(value.(εᴸ)) + sum(value.(εᵁ))

        wᴸ_value = value.(wᴸ); wᵁ_value = value.(wᵁ)
        wᴸ⁻_value = value.(wᴸ⁻); wᵁ⁻_value = value.(wᵁ⁻)
        wᴸ⁺_value = value.(wᴸ⁺); wᵁ⁺_value = value.(wᵁ⁺)
        for i = 1:n
            if wᴸ_value[i] > wᵁ_value[i]
                wᴸ_value[i] = wᵁ_value[i]
            end
            if wᴸ⁻_value[i] > wᵁ⁻_value[i]
                wᴸ⁻_value[i] = wᵁ⁻_value[i]
            end
            if wᴸ⁺_value[i] > wᵁ⁺_value[i]
                wᴸ⁺_value[i] = wᵁ⁺_value[i]
            end
        end
        εᴸ_value = value.(εᴸ); εᵁ_value = value.(εᵁ)

        return (
            wᴸ=wᴸ_value, wᵁ=wᵁ_value,
            wᴸ⁻=wᴸ⁻_value, wᵁ⁻=wᵁ⁻_value,
            wᴸ⁺=wᴸ⁺_value, wᵁ⁺=wᵁ⁺_value,
            εᴸ=εᴸ_value, εᵁ=εᵁ_value,
            optimalValue=optimalValue
        )
    finally
        # エラー終了時にも変数などを消去する
        empty!(model)
    end
end
