using LaTeXStrings

include("./method4.jl")

displayApproximationLPResults_m4 = displayApproximationLPResults_m3

CancellingLPResultLaTeXString_m4 = @NamedTuple{
    W::String, W⁻::String, W⁺::String,
    εᴸ::String, εᵁ::String
    }

function cancellingLPResultLaTeXString_m4(
        result::CancellingLPResult_m4{T}
        )::CancellingLPResultLaTeXString_m4 where {T <: Real}
    n = length(result.wᴸ⁻)

    W = "\\begin{bmatrix}"
    W⁻ = "\\begin{bmatrix}"; W⁺ = "\\begin{bmatrix}"
    εᴸ = "\\begin{bmatrix}"; εᵁ = "\\begin{bmatrix}"
    for i = 1:n
        wᴸ = string(round(result.wᴸ[i], digits=3))
        wᵁ = string(round(result.wᵁ[i], digits=3))
        W *= "\\left[ $(wᴸ), $(wᵁ) \\right]"
        wᴸ⁻ = string(round(result.wᴸ⁻[i], digits=3))
        wᵁ⁻ = string(round(result.wᵁ⁻[i], digits=3))
        W⁻ *= "\\left[ $(wᴸ⁻), $(wᵁ⁻) \\right]"
        wᴸ⁺ = string(round(result.wᴸ⁺[i], digits=3))
        wᵁ⁺ = string(round(result.wᵁ⁺[i], digits=3))
        W⁺ *= "\\left[ $(wᴸ⁺), $(wᵁ⁺) \\right]"
        εᴸ *= " $(string(round(result.εᴸ[i], digits=3))) "
        εᵁ *= " $(string(round(result.εᵁ[i], digits=3))) "
        if i != n
            W *= "\\\\"
            W⁻ *= "\\\\"; W⁺ *= "\\\\"
            εᴸ *= "\\\\"; εᵁ *= "\\\\"
        end
    end
    W *= "\\end{bmatrix}"
    W⁻ *= "\\end{bmatrix}"; W⁺ *= "\\end{bmatrix}"
    εᴸ *= "\\end{bmatrix}"; εᵁ *= "\\end{bmatrix}"

    return (
        W = W, W⁻=W⁻, W⁺=W⁺,
        εᴸ=εᴸ, εᵁ=εᵁ
    )
end

function displayCancellingLPResultLaTeXString_m4(
        result::CancellingLPResult_m4{T}
        ) where {T <: Real}
    resultₛₜᵣ = cancellingLPResultLaTeXString_m4(result)

    display(L"W = %$(resultₛₜᵣ.W)")
    display(L"W^- = %$(resultₛₜᵣ.W⁻)")
    display(L"W^+ = %$(resultₛₜᵣ.W⁺)")
    display(L"""
        \varepsilon^\text{L} = %$(resultₛₜᵣ.εᴸ), ~~
        \varepsilon^\text{U} = %$(resultₛₜᵣ.εᵁ)
        """)
end
