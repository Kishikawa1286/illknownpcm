using JuMP
import HiGHS

include("twofoldIntervalMatrix.jl")
include("utils.jl")

ε = 1e-10

Result_LP_CBA = Tuple{
    # optimal solution
    Vector{T}, # wᵢᵁ
    Vector{T}, # wᵢᴸ
    Vector{T}, # wᵢᴸ⁺
    Vector{T}, # wᵢᴸ⁻
    Vector{T}, # εᵢᴸ
    Vector{T}, # wᵢᵁ⁻
    Vector{T}, # wᵢᵁ⁺
    Vector{T}, # εᵢᵁ
    # optimal value
    T
    } where {T <: Real}

# Twofold Interval PCMが引数
function solveLP_CBA(A::Matrix{TwofoldInterval{T}})::Result_LP_CBA{T} where {T <: Real}
    if !isTwofoldIntervalPCM(A)
        throw(ArgumentError("""
            Given matrix is not valid as twofold interval matrix.
            """))
    end

    m, n = size(A)
    model = Model(HiGHS.Optimizer)
    try
        # 変数の定義
        @variable(model, wᵁ[i=1:n] ≥ ε) # wᵢᵁ ≥ ε
        @variable(model, wᴸ[i=1:n] ≥ ε) # wᵢᴸ ≥ ε
        @variable(model, wᴸ⁺[i=1:n] ≥ 0) # wᵢᴸ⁺ ≥ 0
        @variable(model, wᴸ⁻[i=1:n] ≥ 0) # wᵢᴸ⁻ ≥ ε
        @variable(model, εᴸ[i=1:n] ≥ 0) # εᵢᴸ ≥ 0
        @variable(model, wᵁ⁻[i=1:n] ≥ 0) # wᵢᵁ⁻ ≥ 0
        @variable(model, wᵁ⁺[i=1:n] ≥ 0) # wᵢᵁ⁺ ≥ 0
        @variable(model, εᵁ[i=1:n] ≥ 0) # εᵢᵁ ≥ 0

        for i = 1:n
            wᵢᵁ = wᵁ[i]
            wᵢᴸ = wᴸ[i]
            wᵢᴸ⁺ = wᴸ⁺[i]
            wᵢᴸ⁻ = wᴸ⁻[i]
            εᵢᴸ = εᴸ[i]
            wᵢᵁ⁻ = wᵁ⁻[i]
            wᵢᵁ⁺ = wᵁ⁺[i]
            εᵢᵁ = εᵁ[i]
            
            @constraint(model, wᵢᵁ ≥ wᵢᴸ)
            @constraint(model, εᵢᴸ ≥ wᵢᴸ⁺ - wᵢᴸ⁻)
            @constraint(model, εᵢᵁ ≥ wᵢᵁ⁻ - wᵢᵁ⁺)
        end

        for i = 1:n
            wᵢᵁ = wᵁ[i]
            wᵢᴸ = wᴸ[i]
            wᵢᴸ⁺ = wᴸ⁺[i]
            wᵢᴸ⁻ = wᴸ⁻[i]
            εᵢᴸ = εᴸ[i]
            wᵢᵁ⁻ = wᵁ⁻[i]
            wᵢᵁ⁺ = wᵁ⁺[i]
            εᵢᵁ = εᵁ[i]
            for j = 1:n
                # i ≠ jでの制約
                if i == j
                    continue
                end

                aᵢⱼᴸ⁺ = A[i,j][2].lo
                aᵢⱼᴸ⁻ = A[i,j][1].lo
                aᵢⱼᵁ⁻ = A[i,j][1].hi
                aᵢⱼᵁ⁺ = A[i,j][2].hi

                wⱼᵁ = wᵁ[j]
                wⱼᴸ = wᴸ[j]
                wⱼᴸ⁺ = wᴸ⁺[j]
                wⱼᴸ⁻ = wᴸ⁻[j]
                εⱼᴸ = εᴸ[j]
                wⱼᵁ⁻ = wᵁ⁻[j]
                wⱼᵁ⁺ = wᵁ⁺[j]
                εⱼᵁ = εᵁ[j]

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

        for j = 1:n
            wⱼᵁ = wᵁ[j]
            wⱼᴸ = wᴸ[j]

            # ∑wᵢᵁ + wⱼᴸ ≥ 1
            @constraint(model, sum(map(i -> wᵁ[i], filter(i -> i != j, 1:n))) + wⱼᴸ ≥ 1)
            # ∑wᵢᴸ + wⱼᵁ ≤ 1
            @constraint(model, sum(map(i -> wᴸ[i], filter(i -> i != j, 1:n))) + wⱼᵁ ≤ 1)
        end

        # 目的関数 ∑(εᵢᴸ + εᵢᵁ)
        @objective(model, Min, sum(i -> εᴸ[i] + εᵁ[i], 1:n))

        optimize!(model)
        optimalValue = sum(i -> value(εᴸ[i]) + value(εᵁ[i]), 1:n)

        # Tuple{Vector, Vector, ...}の形で最適値を返す
        # (wᵢᵁ, wᵢᴸ, wᵢᴸ⁺, wᵢᴸ⁻, εᵢᴸ, wᵢᵁ⁻, wᵢᵁ⁺, εᵢᵁ, optimalValue)の順にTupleに入っている
        return (
            # 浮動小数の桁落ちで wᵢᵁ ≥ wᵢᴸ が満たされない場合に修正する
            map(i -> correctPrecisionLoss(value(wᵁ[i]), value(wᴸ[i])), 1:n),
            value.(wᴸ),
            value.(wᴸ⁺),
            value.(wᴸ⁻),
            value.(εᴸ),
            value.(wᵁ⁻),
            value.(wᵁ⁺),
            value.(εᵁ),
            optimalValue
        )
    finally
        empty!(model)
    end
end

function updatePCM_CBA_1(
        A::Matrix{TwofoldInterval{T}},
        lpResult::Result_LP_CBA{T})::Matrix{TwofoldInterval{T}} where {T <: Real}
    if !isTwofoldIntervalPCM(A)
        throw(ArgumentError("""
            Given matrix is not valid as twofold interval matrix.
            """))
    end

    m, n = size(A)
    Â = deepcopy(A)

    wᴸ⁺ = lpResult[3]
    wᴸ⁻ = lpResult[4]
    εᴸ = lpResult[5]
    wᵁ⁻ = lpResult[6]
    wᵁ⁺ = lpResult[7]
    εᵁ = lpResult[8]

    for i = 1:n, j = 1:n
        # 対角成分は (1..1, 1..1) で固定なので更新不要
        if i == j
            continue
        end

        aᵢⱼᴸ⁺ = A[i,j][2].lo
        aᵢⱼᴸ⁻ = A[i,j][1].lo
        aᵢⱼᵁ⁻ = A[i,j][1].hi
        aᵢⱼᵁ⁺ = A[i,j][2].hi

        # âᵢⱼᴸ⁺ = min(aᵢⱼᴸ⁺, (wᵢᴸ⁺ - εᵢᴸ)/(wⱼᵁ⁺ + εⱼᵁ))
        âᵢⱼᴸ⁺ = min(aᵢⱼᴸ⁺, (wᴸ⁺[i] - εᴸ[i])/(wᵁ⁺[j] + εᵁ[j]))
        # âᵢⱼᴸ⁻ = max(aᵢⱼᴸ⁻, (wᵢᴸ⁻ + εᵢᴸ)/(wⱼᵁ⁻ - εⱼᵁ))
        âᵢⱼᴸ⁻ = max(aᵢⱼᴸ⁻, (wᴸ⁻[i] + εᴸ[i])/(wᵁ⁻[j] - εᵁ[j]))
        # âᵢⱼᵁ⁻ = min(aᵢⱼᵁ⁻, (wᵢᵁ⁻ - εᵢᵁ)/(wⱼᴸ⁻ + εⱼᴸ))
        âᵢⱼᵁ⁻ = min(aᵢⱼᵁ⁻, (wᵁ⁻[i] - εᵁ[i])/(wᴸ⁻[j] + εᴸ[j]))
        # âᵢⱼᵁ⁺ = max(aᵢⱼᴸ⁺, (wᵢᵁ⁺ + εᵢᵁ)/(wⱼᴸ⁺ - εⱼᴸ))
        âᵢⱼᵁ⁺ = max(aᵢⱼᵁ⁺, (wᵁ⁺[i] + εᵁ[i])/(wᴸ⁺[j] - εᴸ[j]))

        # âᵢⱼᵁ⁺ = âᵢⱼᴸ⁺ の場合などに桁落ちで âᵢⱼᵁ⁺ < âᵢⱼᴸ⁺ となることがある
        # âᵢⱼᴸ⁺ と âᵢⱼᵁ⁺ が十分に近い値ならば âᵢⱼᴸ⁺ <- âᵢⱼᵁ⁺
        âᵢⱼᴸ⁺ = correctPrecisionLoss(âᵢⱼᴸ⁺, âᵢⱼᵁ⁺)
        # âᵢⱼᴸ⁻ と âᵢⱼᵁ⁻ が十分に近い値ならば âᵢⱼᴸ⁻ <- âᵢⱼᵁ⁻
        âᵢⱼᴸ⁻ = correctPrecisionLoss(âᵢⱼᴸ⁻, âᵢⱼᵁ⁻)

        # (Âᵢⱼ⁻, Âᵢⱼ⁺)
        if âᵢⱼᴸ⁻ > âᵢⱼᵁ⁻
            Â[i, j] = (emptyinterval(), âᵢⱼᴸ⁺..âᵢⱼᵁ⁺)
        else
            Â[i, j] = (âᵢⱼᴸ⁻..âᵢⱼᵁ⁻, âᵢⱼᴸ⁺..âᵢⱼᵁ⁺)
        end
    end

    return Â
end

function updatePCM_CBA_2(
        A::Matrix{TwofoldInterval{T}},
        lpResult::Result_LP_CBA{T})::Matrix{TwofoldInterval{T}} where {T <: Real}
    if !isTwofoldIntervalPCM(A)
        throw(ArgumentError("""
            Given matrix is not valid as twofold interval matrix.
            """))
    end

    m, n = size(A)
    Â = deepcopy(A)

    wᵁ = lpResult[1]
    wᴸ = lpResult[2]
    wᴸ⁺ = lpResult[3]
    wᴸ⁻ = lpResult[4]
    wᵁ⁻ = lpResult[6]
    wᵁ⁺ = lpResult[7]

    for i = 1:n, j = 1:n
        # 対角成分は (1..1, 1..1) で固定なので更新不要
        if i == j
            continue
        end

        aᵢⱼᴸ⁺ = A[i,j][2].lo
        aᵢⱼᴸ⁻ = A[i,j][1].lo
        aᵢⱼᵁ⁻ = A[i,j][1].hi
        aᵢⱼᵁ⁺ = A[i,j][2].hi

        # âᵢⱼᴸ⁺ = min(aᵢⱼᴸ⁺, wᵢᴸ⁻/wⱼᵁ, wᵢᴸ/wⱼᵁ⁻)
        âᵢⱼᴸ⁺ = min(aᵢⱼᴸ⁺, wᴸ⁻[i]/wᵁ[j], wᴸ[i]/wᵁ⁻[j])
        # âᵢⱼᴸ⁻ = max(aᵢⱼᴸ⁻, wᵢᴸ⁺/wⱼᵁ, wᵢᴸ/wⱼᵁ⁺)
        âᵢⱼᴸ⁻ = max(aᵢⱼᴸ⁻, wᴸ⁺[i]/wᵁ[j], wᴸ[i]/wᵁ⁺[j])
        # âᵢⱼᵁ⁻ = min(aᵢⱼᵁ⁻, wᵢᵁ⁺/wⱼᴸ, wᵢᵁ/wⱼᴸ⁺)
        âᵢⱼᵁ⁻ = min(aᵢⱼᵁ⁻, wᵁ⁺[i]/wᴸ[j], wᵁ[i]/wᴸ⁺[j])
        # âᵢⱼᵁ⁺ = max(aᵢⱼᵁ⁺, wᵢᵁ⁻/wⱼᴸ, wᵢᵁ/wⱼᴸ⁻)
        âᵢⱼᵁ⁺ = max(aᵢⱼᵁ⁺, wᵁ⁻[i]/wᴸ[j], wᵁ[i]/wᴸ⁻[j])

        # âᵢⱼᵁ⁺ = âᵢⱼᴸ⁺ の場合などに桁落ちで âᵢⱼᵁ⁺ < âᵢⱼᴸ⁺ となることがある
        # âᵢⱼᴸ⁺ と âᵢⱼᵁ⁺ が十分に近い値ならば âᵢⱼᴸ⁺ <- âᵢⱼᵁ⁺
        âᵢⱼᴸ⁺ = correctPrecisionLoss(âᵢⱼᴸ⁺, âᵢⱼᵁ⁺)
        # âᵢⱼᴸ⁻ と âᵢⱼᵁ⁻ が十分に近い値ならば âᵢⱼᴸ⁻ <- âᵢⱼᵁ⁻
        âᵢⱼᴸ⁻ = correctPrecisionLoss(âᵢⱼᴸ⁻, âᵢⱼᵁ⁻)

        # (Âᵢⱼ⁻, Âᵢⱼ⁺)
        if âᵢⱼᴸ⁻ > âᵢⱼᵁ⁻
            Â[i, j] = (emptyinterval(), âᵢⱼᴸ⁺..âᵢⱼᵁ⁺)
        else
            Â[i, j] = (âᵢⱼᴸ⁻..âᵢⱼᵁ⁻, âᵢⱼᴸ⁺..âᵢⱼᵁ⁺)
        end
    end

    return Â
end
