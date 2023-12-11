using JuMP
import HiGHS

include("../../nearlyEqual/index.jl")
include("../../twofoldInterval/index.jl")
include("../../twofoldIntervalPCM/index.jl")

include("../weight.jl")

function method1(A₁::Matrix{Interval{T}}, A₂::Matrix{Interval{T}})::Matrix{TwofoldInterval{T}} where {T <: Real}
    A = intervalPCM2TwofoldIntervalPCM(A₁, A₂)
    w = twofoldIntervalPCM2Weights(A)
    result = solveLP_m1(A, w.ŵᴸ, w.ŵᵁ)
    Â = updatePCM_m1(A, result)
    return Â
end

LPResult_m1 = @NamedTuple{
    wᴸ::Vector{T}, wᵁ::Vector{T},
    wᴸ⁻::Vector{T}, wᵁ⁻::Vector{T},
    wᴸ⁺::Vector{T}, wᵁ⁺::Vector{T},
    εᴸ::Vector{T}, εᵁ::Vector{T},
    optimalValue::T
    } where {T <: Real}

function solveLP_m1(
        A::Matrix{TwofoldInterval{T}},
        Ŵᴸ::Vector{T},
        Ŵᵁ::Vector{T}
        )::LPResult_m1{T} where {T <: Real}
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
        @variable(model, wᴸ⁻[i=1:n] ≥ 0); @variable(model, wᵁ⁻[i=1:n] ≥ 0)
        # wᵢᴸ⁺ ≥ 0, wᵢᵁ⁺ ≥ 0
        @variable(model, wᴸ⁺[i=1:n] ≥ 0); @variable(model, wᵁ⁺[i=1:n] ≥ 0)
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
        @objective(model, Min, sum(Ŵᴸ' * εᴸ) + sum(Ŵᵁ' * εᵁ))

        optimize!(model)

        optimalValue = sum(value.(Ŵᴸ' * εᴸ)) + sum(value.(Ŵᵁ' * εᵁ))

        return (
            # precision error の補正
            # wᵢᴸ と wᵢᵁ が十分に近い値ならば wᵢᴸ <- wᵢᵁ
            wᴸ=map(i -> correctPrecisionLoss(value(wᴸ[i]), value(wᵁ[i])), 1:n),
            wᵁ=value.(wᵁ),
            # precision error の補正
            # wᵢᴸ⁻ と wᵢᵁ⁻ が十分に近い値ならば wᵢᴸ⁻ <- wᵢᵁ⁻
            wᴸ⁻=map(i -> correctPrecisionLoss(value(wᴸ⁻[i]), value(wᵁ⁻[i])), 1:n),
            wᵁ⁻=value.(wᵁ⁻),
            wᴸ⁺=value.(wᴸ⁺), wᵁ⁺=value.(wᵁ⁺),
            εᴸ=value.(εᴸ), εᵁ=value.(εᵁ),
            optimalValue=optimalValue
        )
    finally
        # エラー終了時にも変数などを消去する
        empty!(model)
    end
end

function updatePCM_m1(
        A::Matrix{TwofoldInterval{T}},
        result::LPResult_m1{T}
        )::Matrix{TwofoldInterval{T}} where {T <: Real}
    if !isTwofoldIntervalPCM(A)
        throw(ArgumentError("Given matrix is not valid as twofold interval matrix."))
    end

    m, n = size(A)
    Â = deepcopy(A)

    wᴸ = result.wᴸ; wᵁ = result.wᵁ
    wᴸ⁻ = result.wᴸ⁻; wᵁ⁻ = result.wᵁ⁻
    wᴸ⁺ = result.wᴸ⁺; wᵁ⁺ = result.wᵁ⁺

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