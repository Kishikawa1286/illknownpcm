include("../crispPCM/index.jl")
include("../evaluation/index.jl")
include("./method1/index.jl")
include("../method2/index.jl")
include("./method3/index.jl")
include("./method4/index.jl")
include("../intervalPCM/index.jl")
include("../twofoldInterval/index.jl")
include("../twofoldIntervalPCM/index.jl")
include("../utils.jl")

function start(
    n::Integer, 
    numOfCrispPCM::Integer, 
    intervalPCMsPerCrispPCM::Integer
    )
    cases = generateSimulationCases(n, numOfCrispPCM, intervalPCMsPerCrispPCM)
    results = Vector{SimulationResult}(undef, length(cases))
    Threads.@threads for i in eachindex(cases)
        results[i] = runSimulation(cases[i])
    end
    conincidenceList = calculateCoincidenceList.(results)
    return total(conincidenceList)
end

# Generate a simulation case

SimulationCase = @NamedTuple{
    A₁::Matrix{Interval{T}},
    A₂::Matrix{Interval{T}}
} where {T <: Real}

function generateSimulationCases(
    n::Integer,
    numOfCrispPCM::Integer,
    intervalPCMsPerCrispPCM::Integer,
    S::Real = 9.0,
    threshold::Real = log(5) / 2
)::Vector{SimulationCase}
    simulation_cases = SimulationCase[]

    for _ in 1:numOfCrispPCM
        # Generate a Crisp PCM
        crisp_pcm = generateConsistentCrispPCM(n, S)
        for _ in 1:intervalPCMsPerCrispPCM
            # Generate an Interval PCM for each Crisp PCM
            interval_pcm_1 = discretizateIntoComparisonScale(randamizedIntervalPCM(crisp_pcm, threshold), S)
            interval_pcm_2 = discretizateIntoComparisonScale(randamizedIntervalPCM(crisp_pcm, threshold), S)
            # Append the Crisp PCM and its corresponding Interval PCM as a tuple to the list
            push!(simulation_cases, (A₁=interval_pcm_1, A₂=interval_pcm_2))
        end
    end

    return simulation_cases
end

# Run each method

SimulationResult = @NamedTuple{
    # 各 DM の与える区間 PCM
    A₁::Matrix{Interval{T}}, A₂::Matrix{Interval{T}},
    # 各 DM の２重区間 PCM
    Ã₁::Matrix{TwofoldInterval{T}}, Ã₂::Matrix{TwofoldInterval{T}},
    Ã₁⁻::Matrix{Interval{T}}, Ã₂⁻::Matrix{Interval{T}},
    Ã₁⁺::Matrix{Interval{T}}, Ã₂⁺::Matrix{Interval{T}},
    # method 1 ~ 4 の二重区間 PCM
    𝓐¹::Matrix{TwofoldInterval{T}}, 𝓐²::Matrix{TwofoldInterval{T}},
    𝓐³::Matrix{TwofoldInterval{T}}, 𝓐⁴::Matrix{TwofoldInterval{T}},
    # 𝓐ᵏ の二重区間の内側の区間
    𝓐¹⁻::Matrix{Interval{T}}, 𝓐²⁻::Matrix{Interval{T}},
    𝓐³⁻::Matrix{Interval{T}}, 𝓐⁴⁻::Matrix{Interval{T}},
    # 𝓐ᵏ の二重区間の外側の区間
    𝓐¹⁺::Matrix{Interval{T}}, 𝓐²⁺::Matrix{Interval{T}},
    𝓐³⁺::Matrix{Interval{T}}, 𝓐⁴⁺::Matrix{Interval{T}},
    # エラーメッセージ
    error::String
} where {T <: Real}

methodList = [method1, method2, method3, method4]

function runSimulation(
    case::SimulationCase{T}
    )::SimulationResult{T} where {T <: Real}
    Ã₁=Ã(case.A₁); Ã₂=Ã(case.A₂)
    Ã₁⁻ = map(Ã₁ᵢⱼ -> Ã₁ᵢⱼ[1], Ã₁)
    Ã₁⁺ = map(Ã₁ᵢⱼ -> Ã₁ᵢⱼ[2], Ã₁)
    Ã₂⁻ = map(Ã₂ᵢⱼ -> Ã₂ᵢⱼ[1], Ã₂)
    Ã₂⁺ = map(Ã₂ᵢⱼ -> Ã₂ᵢⱼ[2], Ã₂)

    𝓐 = Dict(); 𝓐⁻ = Dict(); 𝓐⁺ = Dict()
    
    try
        for k in 1:4
            𝓐[k] = methodList[k](case.A₁, case.A₂)
            𝓐⁻[k] = map(𝓐ᵢⱼ -> 𝓐ᵢⱼ[1], 𝓐[k])
            𝓐⁺[k] = map(𝓐ᵢⱼ -> 𝓐ᵢⱼ[2], 𝓐[k])
        end

        return (
            A₁=case.A₁, A₂=case.A₂,
            Ã₁=Ã₁, Ã₂=Ã₂,
            Ã₁⁻=Ã₁⁻, Ã₂⁻=Ã₂⁻,
            Ã₁⁺=Ã₁⁺, Ã₂⁺=Ã₂⁺,
            𝓐¹=𝓐[1], 𝓐²=𝓐[2], 𝓐³=𝓐[3], 𝓐⁴=𝓐[4],
            𝓐¹⁻=𝓐⁻[1], 𝓐²⁻=𝓐⁻[2], 𝓐³⁻=𝓐⁻[3], 𝓐⁴⁻=𝓐⁻[4],
            𝓐¹⁺=𝓐⁺[1], 𝓐²⁺=𝓐⁺[2], 𝓐³⁺=𝓐⁺[3], 𝓐⁴⁺=𝓐⁺[4],
            error=""
        )
    catch e
        n = size(case.A₁, 2)
        IPCM = fill(1..1, n, n)
        twofoldIPCM = fill((1..1, 1..1), n, n)
        return (
            A₁=case.A₁, A₂=case.A₂,
            Ã₁=Ã₁, Ã₂=Ã₂,
            Ã₁⁻=Ã₁⁻, Ã₂⁻=Ã₂⁻,
            Ã₁⁺=Ã₁⁺, Ã₂⁺=Ã₂⁺,
            𝓐¹=twofoldIPCM, 𝓐²=twofoldIPCM, 𝓐³=twofoldIPCM, 𝓐⁴=twofoldIPCM,
            𝓐¹⁻=IPCM, 𝓐²⁻=IPCM, 𝓐³⁻=IPCM, 𝓐⁴⁻=IPCM,
            𝓐¹⁺=IPCM, 𝓐²⁺=IPCM, 𝓐³⁺=IPCM, 𝓐⁴⁺=IPCM,
            error=string(e)
        )
    end
end

# Calculate the coincidence

InnerConincidenceList = @NamedTuple{
    𝓐¹⁻::T, 𝓐¹⁺::T,
    𝓐²⁻::T, 𝓐²⁺::T,
    𝓐³⁻::T, 𝓐³⁺::T,
    𝓐⁴⁻::T, 𝓐⁴⁺::T
} where {T <: Real}

ConincidenceList = @NamedTuple{
    A₁::InnerConincidenceList{T},
    A₂::InnerConincidenceList{T},
} where {T <: Real}

CoincidenceListTuple = @NamedTuple{
    A::Vector{ConincidenceList{T}}, # ≈, ⊆
    Ã::Vector{ConincidenceList{T}},
    error::String
} where {T <: Real}

function calculateCoincidenceList(
        result::SimulationResult{T}
        )::CoincidenceListTuple{T} where {T <: Real}
    A₁ = result.A₁; A₂ = result.A₂
    Ã₁⁻ = result.Ã₁⁻; Ã₂⁻ = result.Ã₂⁻
    Ã₁⁺ = result.Ã₁⁺; Ã₂⁺ = result.Ã₂⁺
    𝓐¹⁻ = result.𝓐¹⁻; 𝓐¹⁺ = result.𝓐¹⁺
    𝓐²⁻ = result.𝓐²⁻; 𝓐²⁺ = result.𝓐²⁺
    𝓐³⁻ = result.𝓐³⁻; 𝓐³⁺ = result.𝓐³⁺
    𝓐⁴⁻ = result.𝓐⁴⁻; 𝓐⁴⁺ = result.𝓐⁴⁺
    try
        return (
            A=[
                (
                    A₁=(
                        𝓐¹⁻=coincidenceIndex(𝓐¹⁻, A₁), 𝓐¹⁺=coincidenceIndex(A₁, 𝓐¹⁺),
                        𝓐²⁻=coincidenceIndex(𝓐²⁻, A₁), 𝓐²⁺=coincidenceIndex(A₁, 𝓐²⁺),
                        𝓐³⁻=coincidenceIndex(𝓐³⁻, A₁), 𝓐³⁺=coincidenceIndex(A₁, 𝓐³⁺),
                        𝓐⁴⁻=coincidenceIndex(𝓐⁴⁻, A₁), 𝓐⁴⁺=coincidenceIndex(A₁, 𝓐⁴⁺)
                    ),
                    A₂=(
                        𝓐¹⁻=coincidenceIndex(𝓐¹⁻, A₂), 𝓐¹⁺=coincidenceIndex(A₂, 𝓐¹⁺),
                        𝓐²⁻=coincidenceIndex(𝓐²⁻, A₂), 𝓐²⁺=coincidenceIndex(A₂, 𝓐²⁺),
                        𝓐³⁻=coincidenceIndex(𝓐³⁻, A₂), 𝓐³⁺=coincidenceIndex(A₂, 𝓐³⁺),
                        𝓐⁴⁻=coincidenceIndex(𝓐⁴⁻, A₂), 𝓐⁴⁺=coincidenceIndex(A₂, 𝓐⁴⁺)
                    )
                ),
                (
                    A₁=(
                        𝓐¹⁻=containmentIndices(𝓐¹⁻, A₁), 𝓐¹⁺=containmentIndices(A₁, 𝓐¹⁺),
                        𝓐²⁻=containmentIndices(𝓐²⁻, A₁), 𝓐²⁺=containmentIndices(A₁, 𝓐²⁺),
                        𝓐³⁻=containmentIndices(𝓐³⁻, A₁), 𝓐³⁺=containmentIndices(A₁, 𝓐³⁺),
                        𝓐⁴⁻=containmentIndices(𝓐⁴⁻, A₁), 𝓐⁴⁺=containmentIndices(A₁, 𝓐⁴⁺)
                    ),
                    A₂=(
                        𝓐¹⁻=containmentIndices(𝓐¹⁻, A₂), 𝓐¹⁺=containmentIndices(A₂, 𝓐¹⁺),
                        𝓐²⁻=containmentIndices(𝓐²⁻, A₂), 𝓐²⁺=containmentIndices(A₂, 𝓐²⁺),
                        𝓐³⁻=containmentIndices(𝓐³⁻, A₂), 𝓐³⁺=containmentIndices(A₂, 𝓐³⁺),
                        𝓐⁴⁻=containmentIndices(𝓐⁴⁻, A₂), 𝓐⁴⁺=containmentIndices(A₂, 𝓐⁴⁺)
                    )
                )
            ],
            Ã=[
                (
                    A₁=(
                        𝓐¹⁻=coincidenceIndex(𝓐¹⁻, Ã₁⁻), 𝓐¹⁺=coincidenceIndex(Ã₁⁺, 𝓐¹⁺),
                        𝓐²⁻=coincidenceIndex(𝓐²⁻, Ã₁⁻), 𝓐²⁺=coincidenceIndex(Ã₁⁺, 𝓐²⁺),
                        𝓐³⁻=coincidenceIndex(𝓐³⁻, Ã₁⁻), 𝓐³⁺=coincidenceIndex(Ã₁⁺, 𝓐³⁺),
                        𝓐⁴⁻=coincidenceIndex(𝓐⁴⁻, Ã₁⁻), 𝓐⁴⁺=coincidenceIndex(Ã₁⁺, 𝓐⁴⁺)
                    ),
                    A₂=(
                        𝓐¹⁻=coincidenceIndex(𝓐¹⁻, Ã₂⁻), 𝓐¹⁺=coincidenceIndex(Ã₂⁺, 𝓐¹⁺),
                        𝓐²⁻=coincidenceIndex(𝓐²⁻, Ã₂⁻), 𝓐²⁺=coincidenceIndex(Ã₂⁺, 𝓐²⁺),
                        𝓐³⁻=coincidenceIndex(𝓐³⁻, Ã₂⁻), 𝓐³⁺=coincidenceIndex(Ã₂⁺, 𝓐³⁺),
                        𝓐⁴⁻=coincidenceIndex(𝓐⁴⁻, Ã₂⁻), 𝓐⁴⁺=coincidenceIndex(Ã₂⁺, 𝓐⁴⁺)
                    )
                ),
                (
                    A₁=(
                        𝓐¹⁻=containmentIndices(Ã₁⁻, 𝓐¹⁻), 𝓐¹⁺=containmentIndices(𝓐¹⁺, Ã₁⁺),
                        𝓐²⁻=containmentIndices(Ã₁⁻, 𝓐²⁻), 𝓐²⁺=containmentIndices(𝓐²⁺, Ã₁⁺),
                        𝓐³⁻=containmentIndices(Ã₁⁻, 𝓐³⁻), 𝓐³⁺=containmentIndices(𝓐³⁺, Ã₁⁺),
                        𝓐⁴⁻=containmentIndices(Ã₁⁻, 𝓐⁴⁻), 𝓐⁴⁺=containmentIndices(𝓐⁴⁺, Ã₁⁺)
                    ),
                    A₂=(
                        𝓐¹⁻=containmentIndices(Ã₂⁻, 𝓐¹⁻), 𝓐¹⁺=containmentIndices(𝓐¹⁺, Ã₂⁺),
                        𝓐²⁻=containmentIndices(Ã₂⁻, 𝓐²⁻), 𝓐²⁺=containmentIndices(𝓐²⁺, Ã₂⁺),
                        𝓐³⁻=containmentIndices(Ã₂⁻, 𝓐³⁻), 𝓐³⁺=containmentIndices(𝓐³⁺, Ã₂⁺),
                        𝓐⁴⁻=containmentIndices(Ã₂⁻, 𝓐⁴⁻), 𝓐⁴⁺=containmentIndices(𝓐⁴⁺, Ã₂⁺)
                    )
                )
            ],
            error=""
        )
    catch e
        A = (
            𝓐¹⁻=NaN, 𝓐¹⁺=NaN,
            𝓐²⁻=NaN, 𝓐²⁺=NaN,
            𝓐³⁻=NaN, 𝓐³⁺=NaN,
            𝓐⁴⁻=NaN, 𝓐⁴⁺=NaN
        )
        return (
            A=[(A₁=A, A₂=A), (A₁=A, A₂=A)],
            Ã=[(A₁=A, A₂=A), (A₁=A, A₂=A)],
            error=string(e)
        )
    end
end

# Calculate the total coincidence

function total(conincidenceLists::Vector{CoincidenceListTuple{T}}) where {T <: Real}
    if length(conincidenceLists) == 0
        throw(ArgumentError("Empty list of ConincidenceLists"))
    end

    list_A_σ¹⁻ = []; list_A_σ¹⁺ = []
    list_A_σ²⁻ = []; list_A_σ²⁺ = []
    list_A_σ³⁻ = []; list_A_σ³⁺ = []
    list_A_σ⁴⁻ = []; list_A_σ⁴⁺ = []

    list_A_ρ¹⁻ = []; list_A_ρ¹⁺ = []
    list_A_ρ²⁻ = []; list_A_ρ²⁺ = []
    list_A_ρ³⁻ = []; list_A_ρ³⁺ = []
    list_A_ρ⁴⁻ = []; list_A_ρ⁴⁺ = []

    list_Ã_σ¹⁻ = []; list_Ã_σ¹⁺ = []
    list_Ã_σ²⁻ = []; list_Ã_σ²⁺ = []
    list_Ã_σ³⁻ = []; list_Ã_σ³⁺ = []
    list_Ã_σ⁴⁻ = []; list_Ã_σ⁴⁺ = []

    list_Ã_ρ¹⁻ = []; list_Ã_ρ¹⁺ = []
    list_Ã_ρ²⁻ = []; list_Ã_ρ²⁺ = []
    list_Ã_ρ³⁻ = []; list_Ã_ρ³⁺ = []
    list_Ã_ρ⁴⁻ = []; list_Ã_ρ⁴⁺ = []

    for t in conincidenceLists
        if t.error != "" continue end

        push!(list_A_σ¹⁻, sqrt(t.A[1].A₁.𝓐¹⁻ * t.A[1].A₂.𝓐¹⁻))
        push!(list_A_σ¹⁺, sqrt(t.A[1].A₁.𝓐¹⁺ * t.A[1].A₂.𝓐¹⁺))
        push!(list_A_σ²⁻, sqrt(t.A[1].A₁.𝓐²⁻ * t.A[1].A₂.𝓐²⁻))
        push!(list_A_σ²⁺, sqrt(t.A[1].A₁.𝓐²⁺ * t.A[1].A₂.𝓐²⁺))
        push!(list_A_σ³⁻, sqrt(t.A[1].A₁.𝓐³⁻ * t.A[1].A₂.𝓐³⁻))
        push!(list_A_σ³⁺, sqrt(t.A[1].A₁.𝓐³⁺ * t.A[1].A₂.𝓐³⁺))
        push!(list_A_σ⁴⁻, sqrt(t.A[1].A₁.𝓐⁴⁻ * t.A[1].A₂.𝓐⁴⁻))
        push!(list_A_σ⁴⁺, sqrt(t.A[1].A₁.𝓐⁴⁺ * t.A[1].A₂.𝓐⁴⁺))

        push!(list_A_ρ¹⁻, sqrt(t.A[2].A₁.𝓐¹⁻ * t.A[2].A₂.𝓐¹⁻))
        push!(list_A_ρ¹⁺, sqrt(t.A[2].A₁.𝓐¹⁺ * t.A[2].A₂.𝓐¹⁺))
        push!(list_A_ρ²⁻, sqrt(t.A[2].A₁.𝓐²⁻ * t.A[2].A₂.𝓐²⁻))
        push!(list_A_ρ²⁺, sqrt(t.A[2].A₁.𝓐²⁺ * t.A[2].A₂.𝓐²⁺))
        push!(list_A_ρ³⁻, sqrt(t.A[2].A₁.𝓐³⁻ * t.A[2].A₂.𝓐³⁻))
        push!(list_A_ρ³⁺, sqrt(t.A[2].A₁.𝓐³⁺ * t.A[2].A₂.𝓐³⁺))
        push!(list_A_ρ⁴⁻, sqrt(t.A[2].A₁.𝓐⁴⁻ * t.A[2].A₂.𝓐⁴⁻))
        push!(list_A_ρ⁴⁺, sqrt(t.A[2].A₁.𝓐⁴⁺ * t.A[2].A₂.𝓐⁴⁺))

        push!(list_Ã_σ¹⁻, sqrt(t.Ã[1].A₁.𝓐¹⁻ * t.Ã[1].A₂.𝓐¹⁻))
        push!(list_Ã_σ¹⁺, sqrt(t.Ã[1].A₁.𝓐¹⁺ * t.Ã[1].A₂.𝓐¹⁺))
        push!(list_Ã_σ²⁻, sqrt(t.Ã[1].A₁.𝓐²⁻ * t.Ã[1].A₂.𝓐²⁻))
        push!(list_Ã_σ²⁺, sqrt(t.Ã[1].A₁.𝓐²⁺ * t.Ã[1].A₂.𝓐²⁺))
        push!(list_Ã_σ³⁻, sqrt(t.Ã[1].A₁.𝓐³⁻ * t.Ã[1].A₂.𝓐³⁻))
        push!(list_Ã_σ³⁺, sqrt(t.Ã[1].A₁.𝓐³⁺ * t.Ã[1].A₂.𝓐³⁺))
        push!(list_Ã_σ⁴⁻, sqrt(t.Ã[1].A₁.𝓐⁴⁻ * t.Ã[1].A₂.𝓐⁴⁻))
        push!(list_Ã_σ⁴⁺, sqrt(t.Ã[1].A₁.𝓐⁴⁺ * t.Ã[1].A₂.𝓐⁴⁺))

        push!(list_Ã_ρ¹⁻, sqrt(t.Ã[2].A₁.𝓐¹⁻ * t.Ã[2].A₂.𝓐¹⁻))
        push!(list_Ã_ρ¹⁺, sqrt(t.Ã[2].A₁.𝓐¹⁺ * t.Ã[2].A₂.𝓐¹⁺))
        push!(list_Ã_ρ²⁻, sqrt(t.Ã[2].A₁.𝓐²⁻ * t.Ã[2].A₂.𝓐²⁻))
        push!(list_Ã_ρ²⁺, sqrt(t.Ã[2].A₁.𝓐²⁺ * t.Ã[2].A₂.𝓐²⁺))
        push!(list_Ã_ρ³⁻, sqrt(t.Ã[2].A₁.𝓐³⁻ * t.Ã[2].A₂.𝓐³⁻))
        push!(list_Ã_ρ³⁺, sqrt(t.Ã[2].A₁.𝓐³⁺ * t.Ã[2].A₂.𝓐³⁺))
        push!(list_Ã_ρ⁴⁻, sqrt(t.Ã[2].A₁.𝓐⁴⁻ * t.Ã[2].A₂.𝓐⁴⁻))
        push!(list_Ã_ρ⁴⁺, sqrt(t.Ã[2].A₁.𝓐⁴⁺ * t.Ã[2].A₂.𝓐⁴⁺))
    end

    return (
        list_A_σ¹⁻=list_A_σ¹⁻,
        list_A_σ¹⁺=list_A_σ¹⁺,
        list_A_σ²⁻=list_A_σ²⁻,
        list_A_σ²⁺=list_A_σ²⁺,
        list_A_σ³⁻=list_A_σ³⁻,
        list_A_σ³⁺=list_A_σ³⁺,
        list_A_σ⁴⁻=list_A_σ⁴⁻,
        list_A_σ⁴⁺=list_A_σ⁴⁺,
        
        list_A_ρ¹⁻=list_A_ρ¹⁻,
        list_A_ρ¹⁺=list_A_ρ¹⁺,
        list_A_ρ²⁻=list_A_ρ²⁻,
        list_A_ρ²⁺=list_A_ρ²⁺,
        list_A_ρ³⁻=list_A_ρ³⁻,
        list_A_ρ³⁺=list_A_ρ³⁺,
        list_A_ρ⁴⁻=list_A_ρ⁴⁻,
        list_A_ρ⁴⁺=list_A_ρ⁴⁺,
        
        list_Ã_σ¹⁻=list_Ã_σ¹⁻,
        list_Ã_σ¹⁺=list_Ã_σ¹⁺,
        list_Ã_σ²⁻=list_Ã_σ²⁻,
        list_Ã_σ²⁺=list_Ã_σ²⁺,
        list_Ã_σ³⁻=list_Ã_σ³⁻,
        list_Ã_σ³⁺=list_Ã_σ³⁺,
        list_Ã_σ⁴⁻=list_Ã_σ⁴⁻,
        list_Ã_σ⁴⁺=list_Ã_σ⁴⁺,
        
        list_Ã_ρ¹⁻=list_Ã_ρ¹⁻,
        list_Ã_ρ¹⁺=list_Ã_ρ¹⁺,
        list_Ã_ρ²⁻=list_Ã_ρ²⁻,
        list_Ã_ρ²⁺=list_Ã_ρ²⁺,
        list_Ã_ρ³⁻=list_Ã_ρ³⁻,
        list_Ã_ρ³⁺=list_Ã_ρ³⁺,
        list_Ã_ρ⁴⁻=list_Ã_ρ⁴⁻,
        list_Ã_ρ⁴⁺=list_Ã_ρ⁴⁺
    )
end
