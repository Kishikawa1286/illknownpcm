include("../twofoldIntervalPCM/index.jl")

function twofoldIntervalPCM2CrispPCM(
        A::Matrix{TwofoldInterval{T}}
        ) where {T <: Real}
    if !isTwofoldIntervalPCM(A)
        throw(ArgumentError("Given matrix is not valid as twofold interval matrix."))
    end
    
    m, n = size(A)

    Āᴸ = fill(1.0, (n, n)); Āᵁ = fill(1.0, (n, n))
    Āᶜ = fill(1.0, (n, n))

    for i = 1:n, j = 1:n
        if i == j continue end

        aᵢⱼᴸ⁻ = A[i,j][1].lo; aᵢⱼᵁ⁻ = A[i,j][1].hi
        aᵢⱼᴸ⁺ = A[i,j][2].lo; aᵢⱼᵁ⁺ = A[i,j][2].hi

        Āᴸ[i,j] = sqrt(aᵢⱼᴸ⁻ * aᵢⱼᴸ⁺)
        Āᵁ[i,j] = sqrt(aᵢⱼᵁ⁻ * aᵢⱼᵁ⁺)
        Āᶜ[i,j] = sqrt(Āᴸ[i,j] * Āᵁ[i,j])
    end

    return (Āᴸ=Āᴸ, Āᵁ=Āᵁ, Āᶜ=Āᶜ)
end

function geometric_mean(arr::Vector{T}) where {T <: Real}
    product = 1.0
    n = length(arr)
    
    for x in arr
        product *= x
    end
    
    return product^(1/n)
end
