using LaTeXStrings

include("./method1.jl")

LPResultLaTeXString_m1 = @NamedTuple{
    W::String,
    wᴸ⁻::String, wᵁ⁻::String,
    wᴸ⁺::String, wᵁ⁺::String,
    εᴸ::String, εᵁ::String
    }

function lPResultLaTeXString_m1(
        result::LPResult_m1{T}
        )::LPResultLaTeXString_m1 where {T <: Real}
    n = length(result.wᴸ)

    W = "\\begin{bmatrix}"
    wᴸ⁻ = "\\begin{bmatrix}"; wᵁ⁻ = "\\begin{bmatrix}"
    wᴸ⁺ = "\\begin{bmatrix}"; wᵁ⁺ = "\\begin{bmatrix}"
    εᴸ = "\\begin{bmatrix}"; εᵁ = "\\begin{bmatrix}"
    for i = 1:n
        wᴸ = string(round(result.wᴸ[i], digits=3))
        wᵁ = string(round(result.wᵁ[i], digits=3))
        W *= "\\left[ $(wᴸ), $(wᵁ) \\right]"
        wᴸ⁻ *= " $(string(round(result.wᴸ⁻[i], digits=3))) "
        wᵁ⁻ *= " $(string(round(result.wᵁ⁻[i], digits=3))) "
        wᴸ⁺ *= " $(string(round(result.wᴸ⁺[i], digits=3))) "
        wᵁ⁺ *= " $(string(round(result.wᵁ⁺[i], digits=3))) "
        εᴸ *= " $(string(round(result.εᴸ[i], digits=3))) "
        εᵁ *= " $(string(round(result.εᵁ[i], digits=3))) "
        if i != n
            W *= " \\\\ "
            wᴸ⁻ *= " \\\\ "; wᵁ⁻ *= " \\\\ "
            wᴸ⁺ *= " \\\\ "; wᵁ⁺ *= " \\\\ "
            εᴸ *= " \\\\ "; εᵁ *= " \\\\ "
        end
    end
    W *= "\\end{bmatrix}"
    wᴸ⁻ *= "\\end{bmatrix}"; wᵁ⁻ *= "\\end{bmatrix}"
    wᴸ⁺ *= "\\end{bmatrix}"; wᵁ⁺ *= "\\end{bmatrix}"
    εᴸ *= "\\end{bmatrix}"; εᵁ *= "\\end{bmatrix}"

    return (
        W=W, wᴸ⁻=wᴸ⁻, wᵁ⁻=wᵁ⁻,
        wᴸ⁺=wᴸ⁺, wᵁ⁺=wᵁ⁺, εᴸ=εᴸ, εᵁ=εᵁ
    )
end

function displayLPResult_m1(
        result::LPResult_m1{T}
        ) where {T <: Real}
    resultₛₜᵣ = lPResultLaTeXString_m1(result)

    display(L"W = %$(resultₛₜᵣ.W)")
    display(L"""
        w^{\text{L}-} = %$(resultₛₜᵣ.wᴸ⁻), ~~
        w^{\text{U}-} = %$(resultₛₜᵣ.wᵁ⁻), ~~
        w^{\text{L}+} = %$(resultₛₜᵣ.wᴸ⁺), ~~
        w^{\text{U}+} = %$(resultₛₜᵣ.wᵁ⁺)
        """)
    display(L"""
        \varepsilon^\text{L} = %$(resultₛₜᵣ.εᴸ), ~~
        \varepsilon^\text{U} = %$(resultₛₜᵣ.εᵁ)
        """)
end
