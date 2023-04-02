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
            # ∑wᵢᵁ + wⱼᴸ ≥ 1
            wⱼᴸ = wᴸ[j]
            ∑wᵢᵁ = sum(map(i -> wᵁ[i], filter(i -> i != j, 1:n)))
            @constraint(model, ∑wᵢᵁ + wⱼᴸ ≥ 1)
            # ∑wᵢᴸ + wⱼᵁ ≤ 1
            wⱼᵁ = wᵁ[j]
            ∑wᵢᴸ = sum(map(i -> wᴸ[i], filter(i -> i != j, 1:n)))
            @constraint(model, ∑wᵢᴸ + wⱼᵁ ≤ 1)
        end

        # 目的関数 ∑(εᵢᴸ + εᵢᵁ)
        @objective(model, Min, sum(i -> εᴸ[i] + εᵁ[i], 1:n))

        optimize!(model)

        # 最適値の計算
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
        # エラー終了時にも変数などを消去する
        empty!(model)
    end
end

# 重み付き版
function solveLP_CBA_Weighted(A::Matrix{TwofoldInterval{T}})::Result_LP_CBA{T} where {T <: Real}
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
            # ∑wᵢᵁ + wⱼᴸ ≥ 1
            wⱼᴸ = wᴸ[j]
            ∑wᵢᵁ = sum(map(i -> wᵁ[i], filter(i -> i != j, 1:n)))
            @constraint(model, ∑wᵢᵁ + wⱼᴸ ≥ 1)
            # ∑wᵢᴸ + wⱼᵁ ≤ 1
            wⱼᵁ = wᵁ[j]
            ∑wᵢᴸ = sum(map(i -> wᴸ[i], filter(i -> i != j, 1:n)))
            @constraint(model, ∑wᵢᴸ + wⱼᵁ ≤ 1)
        end

        # L+とL-, U+とU-の平均
        āᴸ = zeros(n, n)
        āᵁ = zeros(n, n)
        # ill-knownの4端点の平均
        āᶜ = zeros(n, n)
        for i = 1:n, j = 1:n
            aᵢⱼᴸ⁺ = A[i,j][2].lo
            aᵢⱼᴸ⁻ = A[i,j][1].lo
            aᵢⱼᵁ⁻ = A[i,j][1].hi
            aᵢⱼᵁ⁺ = A[i,j][2].hi
            āᴸ[i,j] = sqrt(aᵢⱼᴸ⁺ * aᵢⱼᴸ⁻)
            āᵁ[i,j] = sqrt(aᵢⱼᵁ⁻ * aᵢⱼᵁ⁺)
            āᶜ[i,j] = sqrt(āᴸ[i,j] * āᵁ[i,j])
        end
        # スケール調整で割るときに使う
        #式の中では変数として置かれていない
        åᶜ = 0
        #幾何平均を取って和を取る
        for j = 1:n
            prod = 1
            for k = 1:n
                prod = prod * āᶜ[j,k]
            end
            åᶜ = (åᶜ)^(1/n) + prod
        end
        display(åᶜ)
        # 重みの計算用変数
        ŵᴸ = ones(n)
        ŵᵁ = ones(n)
        #幾何平均を取って　åᶜ で割る
        for i = 1:n
            for j = 1:n
                ŵᴸ[i] = ŵᴸ[i] * āᴸ[i,j]
                ŵᵁ[i] = ŵᵁ[i] * āᵁ[i,j]
            end
            ŵᴸ[i] = (ŵᴸ[i])^(1/n) / åᶜ
            ŵᵁ[i] = (ŵᵁ[i])^(1/n) / åᶜ
        end

        # 重み（"\{hat\{hat{w}}}"だが文字がないのでẘで代用）
        ẘᴸ = zeros(n)
        ẘᵁ = zeros(n)
        for i = 1:n
            ẘᴸ[i] = max(ŵᴸ[i], 1 - sum(map(j -> ŵᵁ[j], filter(j -> j != i, i:n))))
            ẘᵁ[i] = min(ŵᵁ[i], 1 - sum(map(j -> ŵᴸ[j], filter(j -> j != i, i:n))))
        end

        display(ẘᴸ)
        display(ẘᵁ)
        
        @objective(model, Min, sum(i -> (εᴸ[i] / ẘᴸ[i]) + (εᵁ[i] / ẘᵁ[i]), 1:n))

        optimize!(model)

        # 最適値の計算
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
        # エラー終了時にも変数などを消去する
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

        wᵢᴸ⁺ = wᴸ⁺[i]
        wᵢᴸ⁻ = wᴸ⁻[i]
        wᵢᵁ⁻ = wᵁ⁻[i]
        wᵢᵁ⁺ = wᵁ⁺[i]
        wⱼᴸ⁺ = wᴸ⁺[j]
        wⱼᴸ⁻ = wᴸ⁻[j]
        wⱼᵁ⁻ = wᵁ⁻[j]
        wⱼᵁ⁺ = wᵁ⁺[j]

        εᵢᴸ = εᴸ[i]
        εᵢᵁ = εᵁ[i]
        εⱼᴸ = εᴸ[j]
        εⱼᵁ = εᵁ[j]

        âᵢⱼᴸ⁺ = min(aᵢⱼᴸ⁺, (wᵢᴸ⁺ - εᵢᴸ)/(wⱼᵁ⁺ + εⱼᵁ))
        âᵢⱼᴸ⁻ = max(aᵢⱼᴸ⁻, (wᵢᴸ⁻ + εᵢᴸ)/(wⱼᵁ⁻ - εⱼᵁ))
        âᵢⱼᵁ⁻ = min(aᵢⱼᵁ⁻, (wᵢᵁ⁻ - εᵢᵁ)/(wⱼᴸ⁻ + εⱼᴸ))
        âᵢⱼᵁ⁺ = max(aᵢⱼᵁ⁺, (wᵢᵁ⁺ + εᵢᵁ)/(wⱼᴸ⁺ - εⱼᴸ))

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

        wᵢᵁ = wᵁ[i]
        wᵢᴸ = wᴸ[i]
        wⱼᵁ = wᵁ[j]
        wⱼᴸ = wᴸ[j]

        wᵢᴸ⁺ = wᴸ⁺[i]
        wᵢᴸ⁻ = wᴸ⁻[i]
        wᵢᵁ⁻ = wᵁ⁻[i]
        wᵢᵁ⁺ = wᵁ⁺[i]
        wⱼᴸ⁺ = wᴸ⁺[j]
        wⱼᴸ⁻ = wᴸ⁻[j]
        wⱼᵁ⁻ = wᵁ⁻[j]
        wⱼᵁ⁺ = wᵁ⁺[j]

        âᵢⱼᴸ⁺ = min(aᵢⱼᴸ⁺, wᵢᴸ⁻/wⱼᵁ, wᵢᴸ/wⱼᵁ⁻)
        âᵢⱼᴸ⁻ = max(aᵢⱼᴸ⁻, wᵢᴸ⁺/wⱼᵁ, wᵢᴸ/wⱼᵁ⁺)
        âᵢⱼᵁ⁻ = min(aᵢⱼᵁ⁻, wᵢᵁ⁺/wⱼᴸ, wᵢᵁ/wⱼᴸ⁺)
        âᵢⱼᵁ⁺ = max(aᵢⱼᵁ⁺, wᵢᵁ⁻/wⱼᴸ, wᵢᵁ/wⱼᴸ⁻)

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
