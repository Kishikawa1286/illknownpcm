using IntervalArithmetic
using JuMP
using NLopt
using Ipopt

include("./solution.jl")
include("../twofoldInterval/index.jl")
include("../twofoldIntervalPCM/index.jl")

function solveLogarithmic(
        A::Matrix{TwofoldInterval{T}}
        )::Solution{T} where {T <: Real}
    ε = 1e-8

    if !isTwofoldIntervalPCM(A)
        throw(ArgumentError("Given matrix is not valid as twofold interval matrix."))
    end

    m, n = size(A)

    # 初期値
    w̄ᴸ = fill(1/n, n); w̄ᵁ = fill(1/n, n)
    # Vector の入れ物を準備
    w̄ᴸ⁻ = fill(ε, n); w̄ᵁ⁻ = fill(ε, n)
    w̄ᴸ⁺ = fill(ε, n); w̄ᵁ⁺ = fill(ε, n)
    ε̄¯ᴸ = fill(0.0, n); ε̄¯ᵁ = fill(0.0, n)
    for i = 1:n, j = 1:n
        aᵢⱼᴸ⁻ = A[i,j][1].lo; aᵢⱼᵁ⁻ = A[i,j][1].hi
        aᵢⱼᴸ⁺ = A[i,j][2].lo; aᵢⱼᵁ⁺ = A[i,j][2].hi

        if aᵢⱼᴸ⁻ / n < w̄ᴸ⁻[i]
            w̄ᴸ⁻[i] = aᵢⱼᴸ⁻ / n
        end
        if aᵢⱼᵁ⁻ / n > w̄ᵁ⁻[i]
            w̄ᵁ⁻[i] = aᵢⱼᵁ⁻ / n
        end
        if aᵢⱼᴸ⁺ / n > w̄ᴸ⁺[i]
            w̄ᴸ⁺[i] = aᵢⱼᴸ⁺ / n
        end
        if aᵢⱼᵁ⁺ / n < w̄ᵁ⁺[i]
            w̄ᵁ⁺[i] = aᵢⱼᵁ⁺ / n
        end

        aⱼᵢᵁ⁻ = A[j,i][1].hi; aⱼᵢᴸ⁺ = A[j,i][2].lo

        if log(aᵢⱼᴸ⁺) > ε̄¯ᴸ[i]
            ε̄¯ᴸ[i] = log(aᵢⱼᴸ⁺)
        end
        if log(aⱼᵢᵁ⁻) > ε̄¯ᴸ[i]
            ε̄¯ᴸ[i] = log(aⱼᵢᵁ⁻)
        end
        if log(aᵢⱼᵁ⁻) > ε̄¯ᵁ[i]
            ε̄¯ᵁ[i] = log(aᵢⱼᵁ⁻)
        end
        if log(aⱼᵢᴸ⁺) > ε̄¯ᵁ[i]
            ε̄¯ᵁ[i] = log(aⱼᵢᴸ⁺)
        end
    end

    # model = Model(NLopt.Optimizer)
    # # 参照: https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/
    # set_optimizer_attribute(model, "algorithm", :AUGLAG)
    # set_optimizer_attribute(model, "local_optimizer", :LD_MMA)
    # # set_optimizer_attribute のコメントを参照
    # set_optimizer_attribute(model, "tol", 1e-4)
    # set_optimizer_attribute(model, "constrtol_abs", 1e-4)
    # set_optimizer_attribute(model, "max_iter", 100)

    model = Model(Ipopt.Optimizer)
    set_attribute(model, "max_cpu_time", 60.0)
    set_attribute(model, "print_level", 0)
    
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

            @NLconstraint(model, wᵢᵁ ≥ wᵢᴸ)
            @NLconstraint(model, εᵢᴸ ≥ log(wᵢᴸ⁺) - log(wᵢᴸ⁻))
            @NLconstraint(model, εᵢᵁ ≥ log(wᵢᵁ⁻) - log(wᵢᵁ⁺))

            # 正規性条件
            ∑wⱼᴸ = zero(T); ∑wⱼᵁ = zero(T)
            for j = 1:n
                if i == j
                    continue
                end
                ∑wⱼᴸ += wᴸ[j]; ∑wⱼᵁ += wᵁ[j]
            end
            @NLconstraint(model, ∑wⱼᴸ + wᵢᵁ ≤ 1)
            @NLconstraint(model, ∑wⱼᵁ + wᵢᴸ ≥ 1)

            for j = 1:n
                if i == j
                    continue
                end

                aᵢⱼᴸ⁻ = A[i,j][1].lo; aᵢⱼᵁ⁻ = A[i,j][1].hi
                aᵢⱼᴸ⁺ = A[i,j][2].lo; aᵢⱼᵁ⁺ = A[i,j][2].hi
                wⱼᴸ = wᴸ[j]; wⱼᵁ = wᵁ[j]

                @NLconstraint(model, log(aᵢⱼᴸ⁺) + log(wⱼᵁ) - εᵢᴸ ≤ log(wᵢᴸ))
                @NLconstraint(model, log(wᵢᴸ) ≤ log(aᵢⱼᴸ⁻) + log(wⱼᵁ) + εᵢᴸ)
                @NLconstraint(model, log(aᵢⱼᵁ⁻) + log(wⱼᴸ) - εᵢᵁ ≤ log(wᵢᵁ))
                @NLconstraint(model, log(wᵢᵁ) ≤ log(aᵢⱼᵁ⁺) + log(wⱼᴸ) + εᵢᵁ)

                @NLconstraint(model, wᵢᴸ⁻ ≤ aᵢⱼᴸ⁻ * wⱼᵁ)
                @NLconstraint(model, wᵢᴸ⁺ ≥ aᵢⱼᴸ⁺ * wⱼᵁ)
                @NLconstraint(model, wᵢᵁ⁻ ≥ aᵢⱼᵁ⁻ * wⱼᴸ)
                @NLconstraint(model, wᵢᵁ⁺ ≤ aᵢⱼᵁ⁺ * wⱼᴸ)
            end
        end
        Σwᴸ = zero(T); Σwᵁ = zero(T)
        for i = 1:n
            Σwᴸ += wᴸ[i]; Σwᵁ += wᵁ[i]
        end
        @NLconstraint(model, Σwᴸ + Σwᵁ == 2)
        

        # 目的関数 ∑(εᵢᴸ + εᵢᵁ)
        obj = zero(T)
        for i = 1:n
            obj += εᴸ[i] + εᵁ[i]
        end
        @NLobjective(model, Min, obj)

        for i = 1:n
            set_start_value(wᴸ[i], w̄ᴸ[i]); set_start_value(wᵁ[i], w̄ᵁ[i])
            set_start_value(wᴸ⁻[i], w̄ᴸ⁻[i]); set_start_value(wᵁ⁻[i], w̄ᵁ⁻[i])
            set_start_value(wᴸ⁺[i], w̄ᴸ⁺[i]); set_start_value(wᵁ⁺[i], w̄ᵁ⁺[i])
            set_start_value(εᴸ[i], ε̄¯ᴸ[i]); set_start_value(εᵁ[i], ε̄¯ᵁ[i])
        end

        JuMP.optimize!(model)

        optimalValue = sum(value.(εᴸ)) + sum(value.(εᵁ))

        wᴸ_value = value.(wᴸ); wᵁ_value = value.(wᵁ)
        wᴸ⁻_value = value.(wᴸ⁻); wᵁ⁻_value = value.(wᵁ⁻)
        for i = 1:n
            if wᴸ_value[i] > wᵁ_value[i]
                wᴸ_value[i] = wᵁ_value[i]
            end
            if wᴸ⁻_value[i] > wᵁ⁻_value[i]
                wᴸ⁻_value[i] = wᵁ⁻_value[i]
            end
        end
        wᴸ⁺_value = value.(wᴸ⁺); wᵁ⁺_value = value.(wᵁ⁺)
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
