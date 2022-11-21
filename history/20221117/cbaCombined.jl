using JuMP
import HiGHS

include("twofoldIntervalMatrix.jl")
include("utils.jl")

ε = 1e-10

Result_LP_CBA_Combined = Tuple{
    # optimal solution
    Vector{T}, # wᵢᵁ
    Vector{T}, # wᵢᴸ
    Matrix{T}, # εᵢⱼᴸ⁻
    Matrix{T}, # εᵢⱼᵁ⁻
    Matrix{T}, # εᵢⱼᴸ⁺
    Matrix{T}, # εᵢⱼᵁ⁺
    # optimal value
    T
    } where {T <: Real}

# Twofold Interval PCMが引数
function solveLP_CBA_Combined(
        A::Matrix{TwofoldInterval{T}})::Result_LP_CBA_Combined{T} where {T <: Real}
    if !isTwofoldIntervalPCM(A)
        throw(ArgumentError("""
            Given matrix is not valid as twofold interval matrix.
            """))
    end

    m, n = size(A)
    model = Model(HiGHS.Optimizer)
    try
        # 変数の定義
        @variable(model, wᵁ[i=1:n] >= ε) # wᵢᵁ ≥ ε
        @variable(model, wᴸ[i=1:n] >= ε) # wᵢᴸ ≥ ε
        # Matrixに入れるのでLPには存在しないi==jの変数も定義する
        # 目的関数を最小化する影響でこれらは0になる
        @variable(model, εᴸ⁻[i=1:n,j=1:n] >= 0) # εᵢⱼᴸ⁻ ≥ 0
        @variable(model, εᵁ⁻[i=1:n,j=1:n] >= 0) # εᵢⱼᵁ⁻ ≥ 0
        @variable(model, εᴸ⁺[i=1:n,j=1:n] >= 0) # εᵢⱼᴸ⁺ ≥ 0
        @variable(model, εᵁ⁺[i=1:n,j=1:n] >= 0) # εᵢⱼᵁ⁺ ≥ 0

        for i = 1:n
            wᵢᵁ = wᵁ[i]
            wᵢᴸ = wᴸ[i]
            @constraint(model, wᵢᵁ >= wᵢᴸ)
        end

        for i = 1:n
            wᵢᵁ = wᵁ[i]
            wᵢᴸ = wᴸ[i]
            for j = 1:n
                # i ≠ jでの制約
                if i == j
                    continue
                end

                wⱼᵁ = wᵁ[j]
                wⱼᴸ = wᴸ[j]
                aᵢⱼᴸ⁺ = A[i,j][2].lo
                aᵢⱼᴸ⁻ = A[i,j][1].lo
                aᵢⱼᵁ⁻ = A[i,j][1].hi
                aᵢⱼᵁ⁺ = A[i,j][2].hi
                εᵢⱼᴸ⁻ = εᴸ⁻[i,j]
                εᵢⱼᵁ⁻ = εᵁ⁻[i,j]
                εᵢⱼᴸ⁺ = εᴸ⁺[i,j]
                εᵢⱼᵁ⁺ = εᵁ⁺[i,j]

                @constraint(model, aᵢⱼᴸ⁺ * wⱼᵁ - εᵢⱼᴸ⁺ <= wᵢᴸ)
                @constraint(model, wᵢᴸ <= aᵢⱼᴸ⁻ * wⱼᵁ + εᵢⱼᴸ⁻)
                @constraint(model, aᵢⱼᵁ⁻ * wⱼᴸ - εᵢⱼᵁ⁻ <= wᵢᵁ)
                @constraint(model, wᵢᵁ <= aᵢⱼᵁ⁺ * wⱼᴸ + εᵢⱼᵁ⁺)
            end
        end

        for j = 1:n
            # ∑wᵢᵁ + wⱼᴸ ≥ 1
            @constraint(model, sum(map(i -> wᵁ[i], filter(i -> i != j, 1:n))) + wᴸ[j] >= 1)
            # ∑wᵢᴸ + wⱼᵁ ≤ 1
            @constraint(model, sum(map(i -> wᴸ[i], filter(i -> i != j, 1:n))) + wᵁ[j] <= 1)
        end

        # 目的関数 ∑(εᵢᴸ + εᵢᵁ)
        objective = 0
        for i = 1:n, j = 1:n
            # 値を0にするためにi==jも入れる
            εᵢⱼᴸ⁻ = εᴸ⁻[i,j]
            εᵢⱼᵁ⁻ = εᵁ⁻[i,j]
            εᵢⱼᴸ⁺ = εᴸ⁺[i,j]
            εᵢⱼᵁ⁺ = εᵁ⁺[i,j]
            objective = objective + εᵢⱼᴸ⁻ + εᵢⱼᵁ⁻ + εᵢⱼᴸ⁺ + εᵢⱼᵁ⁺
        end
        @objective(model, Min, objective)

        optimize!(model)

        # 最適値の計算
        optimalValue = 0
        for i = 1:n, j = 1:n
            if i == j
                continue
            end
            εᵢⱼᴸ⁻ = value(εᴸ⁻[i,j])
            εᵢⱼᵁ⁻ = value(εᵁ⁻[i,j])
            εᵢⱼᴸ⁺ = value(εᴸ⁺[i,j])
            εᵢⱼᵁ⁺ = value(εᵁ⁺[i,j])
            optimalValue = optimalValue + εᵢⱼᴸ⁻ + εᵢⱼᵁ⁻ + εᵢⱼᴸ⁺ + εᵢⱼᵁ⁺
        end

        # Tuple{Vector, Vector, ...}の形で最適値を返す
        # (wᵢᵁ, wᵢᴸ, wᵢᴸ⁺, wᵢᴸ⁻, εᵢᴸ, wᵢᵁ⁻, wᵢᵁ⁺, εᵢᵁ, optimalValue)の順にTupleに入っている
        return (
            # 浮動小数の桁落ちで wᵢᵁ ≥ wᵢᴸ が満たされない場合に修正する
            map(i -> correctPrecisionLoss(value(wᵁ[i]), value(wᴸ[i])), 1:n),
            value.(wᴸ),
            value.(εᴸ⁻),
            value.(εᵁ⁻),
            value.(εᴸ⁺),
            value.(εᵁ⁺),
            optimalValue
        )
    finally
        empty!(model)
    end
end

function updatePCM_CBA_Combined(
        A::Matrix{TwofoldInterval{T}},
        lpResult::Result_LP_CBA_Combined{T})::Matrix{TwofoldInterval{T}} where {T <: Real}
    if !isTwofoldIntervalPCM(A)
        throw(ArgumentError("""
            Given matrix is not valid as twofold interval matrix.
            """))
    end

    m, n = size(A)
    Â = deepcopy(A)

    wᵁ = lpResult[1]
    wᴸ = lpResult[2]
    εᴸ⁻ = lpResult[3]
    εᵁ⁻ = lpResult[4]
    εᴸ⁺ = lpResult[5]
    εᵁ⁺ = lpResult[6]

    for i = 1:n
        wᵢᵁ = wᵁ[i]
        wᵢᴸ = wᴸ[i]
        for j = 1:n
            # 対角成分は (1..1, 1..1) で固定なので更新不要
            if i == j
                continue
            end

            wⱼᵁ = wᵁ[j]
            wⱼᴸ = wᴸ[j]

            aᵢⱼᴸ⁺ = A[i,j][2].lo
            aᵢⱼᴸ⁻ = A[i,j][1].lo
            aᵢⱼᵁ⁻ = A[i,j][1].hi
            aᵢⱼᵁ⁺ = A[i,j][2].hi
            aⱼᵢᴸ⁺ = A[j,i][2].lo
            aⱼᵢᴸ⁻ = A[j,i][1].lo
            aⱼᵢᵁ⁻ = A[j,i][1].hi
            aⱼᵢᵁ⁺ = A[j,i][2].hi

            εᵢⱼᴸ⁻ = εᴸ⁻[i,j]
            εᵢⱼᵁ⁻ = εᵁ⁻[i,j]
            εᵢⱼᴸ⁺ = εᴸ⁺[i,j]
            εᵢⱼᵁ⁺ = εᵁ⁺[i,j]
            εⱼᵢᴸ⁻ = εᴸ⁻[j,i]
            εⱼᵢᵁ⁻ = εᵁ⁻[j,i]
            εⱼᵢᴸ⁺ = εᴸ⁺[j,i]
            εⱼᵢᵁ⁺ = εᵁ⁺[j,i]

            âᵢⱼᴸ⁺ = min(aᵢⱼᴸ⁺ - εᵢⱼᴸ⁺/wⱼᵁ, (aⱼᵢᵁ⁺ + εⱼᵢᵁ⁺/wᵢᴸ)^(-1))
            âᵢⱼᴸ⁻ = max(aᵢⱼᴸ⁻ + εᵢⱼᴸ⁻/wⱼᵁ, (aⱼᵢᵁ⁻ - εⱼᵢᵁ⁻/wᵢᴸ)^(-1))
            âᵢⱼᵁ⁻ = min(aᵢⱼᵁ⁻ - εᵢⱼᵁ⁻/wⱼᴸ, (aⱼᵢᴸ⁻ + εⱼᵢᴸ⁻/wᵢᵁ)^(-1))
            âᵢⱼᵁ⁺ = max(aᵢⱼᵁ⁺ + εᵢⱼᵁ⁺/wⱼᴸ, (aⱼᵢᴸ⁺ - εⱼᵢᴸ⁺/wᵢᵁ)^(-1))

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
    end

    return Â
end
