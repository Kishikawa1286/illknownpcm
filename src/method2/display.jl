using LaTeXStrings

include("./method2.jl")

ApproximationLPResultLaTeXString_m2 = @NamedTuple{
    W⁻::String, W⁺::String,
    εᴸ⁻::String, εᵁ⁻::String,
    εᴸ⁺::String, εᵁ⁺::String
    }

function approximationLPResultLaTeXString_m2(
        result::ApproximationLPResult_m2{T}
        )::ApproximationLPResultLaTeXString_m2 where {T <: Real}
    n = length(result.wᴸ⁻)

    W⁻ = "\\begin{bmatrix}"; W⁺ = "\\begin{bmatrix}"
    εᴸ⁻ = "\\begin{bmatrix}"; εᵁ⁻ = "\\begin{bmatrix}"
    εᴸ⁺ = "\\begin{bmatrix}"; εᵁ⁺ = "\\begin{bmatrix}"
    for i = 1:n
        wᴸ⁻ = string(round(result.wᴸ⁻[i], digits=3))
        wᵁ⁻ = string(round(result.wᵁ⁻[i], digits=3))
        W⁻ *= "\\left[ $(wᴸ⁻), $(wᵁ⁻) \\right]"
        wᴸ⁺ = string(round(result.wᴸ⁺[i], digits=3))
        wᵁ⁺ = string(round(result.wᵁ⁺[i], digits=3))
        W⁺ *= "\\left[ $(wᴸ⁺), $(wᵁ⁺) \\right]"
        εᴸ⁻ *= " $(string(round(result.εᴸ⁻[i], digits=3))) "
        εᵁ⁻ *= " $(string(round(result.εᵁ⁻[i], digits=3))) "
        εᴸ⁺ *= " $(string(round(result.εᴸ⁺[i], digits=3))) "
        εᵁ⁺ *= " $(string(round(result.εᵁ⁺[i], digits=3))) "
        if i != n
            W⁻ *= "\\\\"; W⁺ *= "\\\\"
            εᴸ⁻ *= "\\\\"; εᵁ⁻ *= "\\\\"
            εᴸ⁺ *= "\\\\"; εᵁ⁺ *= "\\\\"
        end
    end
    W⁻ *= "\\end{bmatrix}"; W⁺ *= "\\end{bmatrix}"
    εᴸ⁻ *= "\\end{bmatrix}"; εᵁ⁻ *= "\\end{bmatrix}"
    εᴸ⁺ *= "\\end{bmatrix}"; εᵁ⁺ *= "\\end{bmatrix}"

    return (
        W⁻=W⁻, W⁺=W⁺,
        εᴸ⁻=εᴸ⁻, εᵁ⁻=εᵁ⁻, εᴸ⁺=εᴸ⁺, εᵁ⁺=εᵁ⁺
    )
end

function displayApproximationLPResult_m2(
        result::ApproximationLPResult_m2{T}
        ) where {T <: Real}
    resultₛₜᵣ = approximationLPResultLaTeXString_m2(result)

    display(L"W^- = %$(resultₛₜᵣ.W⁻)")
    display(L"W^+ = %$(resultₛₜᵣ.W⁺)")
    display(L"""
        \varepsilon^{\text{L}-} = %$(resultₛₜᵣ.εᴸ⁻), ~~
        \varepsilon^{\text{U}-} = %$(resultₛₜᵣ.εᵁ⁻), ~~
        \varepsilon^{\text{L}+} = %$(resultₛₜᵣ.εᴸ⁺), ~~
        \varepsilon^{\text{U}+} = %$(resultₛₜᵣ.εᵁ⁺)
        """)
end
