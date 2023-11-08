using IntervalArithmetic
using JuMP
import HiGHS

include("../intervalPCM/index.jl")

UpperSolution = @NamedTuple{
    wᴸ::Vector{T}, wᵁ::Vector{T},
    optimalValue::T
    } where {T <: Real}

function solveUpper(
        A::Matrix{Interval{T}}
        )::UpperSolution{T} where {T <: Real}
    ε = 1e-8

    if !isIntervalPCM(A)
        throw(ArgumentError("Given matrix is not valid as interval PCM."))
    end

    m, n = size(A)
    model = Model(HiGHS.Optimizer)
    set_silent(model)

    try
        @variable(model, wᴸ[i=1:n] ≥ ε); @variable(model, wᵁ[i=1:n] ≥ ε)
        
        for i = 1:n
            wᵢᴸ = wᴸ[i]; wᵢᵁ = wᵁ[i]

            @constraint(model, wᵢᵁ ≥ wᵢᴸ)

            ∑wⱼᴸ = sum(map(j -> wᴸ[j], filter(j -> i != j, 1:n)))
            @constraint(model, ∑wⱼᴸ + wᵢᵁ ≤ 1)
            ∑wⱼᵁ = sum(map(j -> wᵁ[j], filter(j -> i != j, 1:n)))
            @constraint(model, ∑wⱼᵁ + wᵢᴸ ≥ 1)

            for j = 1:n
                if i == j continue end

                aᵢⱼᴸ = A[i,j].lo; aᵢⱼᵁ = A[i,j].hi
                wⱼᴸ = wᴸ[j]; wⱼᵁ = wᵁ[j]

                @constraint(model, wᵢᴸ ≤ aᵢⱼᴸ * wⱼᵁ)
                @constraint(model, aᵢⱼᵁ * wⱼᴸ ≤ wᵢᵁ)
            end
        end

        @constraint(model, sum(wᴸ) + sum(wᵁ) == 2)

        @objective(model, Min, sum(wᵁ) - sum(wᴸ))

        optimize!(model)

        optimalValue = sum(value.(wᴸ)) + sum(value.(wᵁ))

        wᴸ_value = value.(wᴸ); wᵁ_value = value.(wᵁ)
        for i = 1:n
            if wᴸ_value[i] > wᵁ_value[i]
                wᴸ_value[i] = wᵁ_value[i]
            end
        end

        return (
            wᴸ = wᴸ_value,
            wᵁ = wᵁ_value,
            optimalValue = optimalValue
        )
    finally
        empty!(model)
    end
end
