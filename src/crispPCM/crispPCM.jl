using LinearAlgebra
using Random
using Distributions

@inline function isCrispPCM(A::Matrix{T})::Bool where {T <: Real}
    # Check if the matrix is square
    if size(A, 1) != size(A, 2)
        return false
    end

    # Check reciprocity
    for i in 1:size(A, 1)
        for j in (i+1):size(A, 2)
            if !nearlyEqual(A[i, j], 1/A[j, i])
                return false
            end
        end
    end

    # Check if the diagonal elements are 1
    if any(diag(A) .!= 1)
        return false
    end

    return true
end

@inline function generateConsistentCrispPCM(
    n::Integer,
    S::T,
    max_tries::Integer = 1000
    )::Matrix{T} where {T <: Real}
    min_cr = Inf
    best_A = generateCrispPCM(n, S) 

    for _ in 1:max_tries
        A = generateCrispPCM(n, S)
        cr = consistencyRatio(A)
        if cr < 0.1
            return A
        elseif cr < min_cr
            min_cr = cr
            best_A = A
        end
    end

    return best_A
end

@inline function generateCrispPCM(n::Integer, S::T)::Matrix where {T <: Real}
    if S <= 1
        throw(ArgumentError("S must be larger than 1."))
    end

    A = ones(n, n)

    for i in 1:n
        for j in (i+1):n
            aᵢⱼ = exp(rand(Uniform(-log(S), log(S))))
            A[i, j] = aᵢⱼ
            A[j, i] = 1/aᵢⱼ
        end
    end

    return A
end

@inline function randomIndex(n::Integer)::Real
    if n < 3
        return 0
    end
    if n == 3
        return 0.58
    end
    if n == 4
        return 0.9
    end
    if n == 5
        return 1.12
    end
    if n == 6
        return 1.24
    end
    if n == 7
        return 1.32
    end
    if n == 8
        return 1.41
    end
    if n == 9
        return 1.45
    end
    if n >= 10
        return 1.49
    end
end

@inline function consistencyIndex(A::Matrix{T})::Real where {T <: Real}
    if (!isCrispPCM(A))
        throw(ArgumentError("A must be crisp PCM."))
    end

    n = size(A, 2)
    λₘₐₓ = maximum(real.(filter(λ  -> isreal(λ), eigvals(A))))
    return (λₘₐₓ - n) / (n - 1)
end


@inline function consistencyRatio(A::Matrix{T})::Real where {T <: Real}
    if (!isCrispPCM(A))
        throw(ArgumentError("A must be crisp PCM."))
    end

    n = size(A, 2)
    return consistencyIndex(A) / randomIndex(n)
end
