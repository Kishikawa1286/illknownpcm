using LaTeXStrings

include("../concatImportanceMethod/concatImportanceMethod.jl")

ConcatImportanceMethodApproximationLPResultLaTeXString = @NamedTuple{
    wₖᴸ⁻::String, wₖᵁ⁻::String,
    wₖᴸ⁺::String, wₖᵁ⁺::String
    }

function concatImportanceMethodApproximationLPResultLaTeXString(
        result::ConcatImportanceMethodApproximationLPResult{T}
        )::ConcatImportanceMethodApproximationLPResultLaTeXString where {T <: Real}
    n = length(result.wₖᴸ⁻)

    wₖᴸ⁻ = "\\begin{bmatrix}"; wₖᵁ⁻ = "\\begin{bmatrix}"
    wₖᴸ⁺ = "\\begin{bmatrix}"; wₖᵁ⁺ = "\\begin{bmatrix}"
    for i = 1:n
        wₖᴸ⁻ *= " $(string(round(result.wₖᴸ⁻[i], digits=3))) "
        wₖᵁ⁻ *= " $(string(round(result.wₖᵁ⁻[i], digits=3))) "
        wₖᴸ⁺ *= " $(string(round(result.wₖᴸ⁺[i], digits=3))) "
        wₖᵁ⁺ *= " $(string(round(result.wₖᵁ⁺[i], digits=3))) "
        if i != n
            wₖᴸ⁻ *= " \\\\ "; wₖᵁ⁻ *= " \\\\ "
            wₖᴸ⁺ *= " \\\\ "; wₖᵁ⁺ *= " \\\\ "
        end
    end
    wₖᴸ⁻ *= "\\end{bmatrix}"; wₖᵁ⁻ *= "\\end{bmatrix}"
    wₖᴸ⁺ *= "\\end{bmatrix}"; wₖᵁ⁺ *= "\\end{bmatrix}"

    return (wₖᴸ⁻=wₖᴸ⁻, wₖᵁ⁻=wₖᵁ⁻, wₖᴸ⁺=wₖᴸ⁺, wₖᵁ⁺=wₖᵁ⁺)
end

function displayConcatImportanceMethodApproximationLPResults(
        results::AbstractArray{ConcatImportanceMethodApproximationLPResult{T}}
        ) where {T <: Real}
    for k = eachindex(results)
        resultₛₜᵣ = concatImportanceMethodApproximationLPResultLaTeXString(results[k])

        display(L"""
            w_{%$(k)}^{\text{L}-} = %$(resultₛₜᵣ.wₖᴸ⁻), ~~
            w_{%$(k)}^{\text{U}-} = %$(resultₛₜᵣ.wₖᵁ⁻), ~~
            w_{%$(k)}^{\text{L}+} = %$(resultₛₜᵣ.wₖᴸ⁺), ~~
            w_{%$(k)}^{\text{U}+} = %$(resultₛₜᵣ.wₖᵁ⁺)
        """)
    end
end

ConcatImportanceMethodTBoundariesLaTeXString = @NamedTuple{
    tᴸ⁻::String, tᵁ⁻::String, tᴸ⁺::String, tᵁ⁺::String
    }

function concatImportanceMethodTBoundariesLaTeXString(
        boundaries::ConcatImportanceMethodTBoundaries{T}
        )::AbstractArray{ConcatImportanceMethodTBoundariesLaTeXString} where {T <: Real}
    m = length(boundaries.tᴸ⁻)

    boundariesₛₜᵣ = fill((tᴸ⁻="", tᵁ⁻="", tᴸ⁺="", tᵁ⁺=""), m)

    for k = 1:m
        tᴸ⁻=""; tᵁ⁻=""; tᴸ⁺=""; tᵁ⁺=""
        tᴸ⁻ *= " $(string(round(boundaries.tᴸ⁻[k], digits=3))) "
        tᵁ⁻ *= " $(string(round(boundaries.tᵁ⁻[k], digits=3))) "
        tᴸ⁺ *= " $(string(round(boundaries.tᴸ⁺[k], digits=3))) "
        tᵁ⁺ *= " $(string(round(boundaries.tᵁ⁺[k], digits=3))) "
        boundariesₛₜᵣ[k] = (tᴸ⁻=tᴸ⁻, tᵁ⁻=tᵁ⁻, tᴸ⁺=tᴸ⁺, tᵁ⁺=tᵁ⁺)
    end

    return boundariesₛₜᵣ
end

function displayConcatImportanceMethodTBoundaries(
        boundaries::ConcatImportanceMethodTBoundaries{T}) where {T <: Real}
    boundariesₛₜᵣ = concatImportanceMethodTBoundariesLaTeXString(boundaries)

    for k = eachindex(boundariesₛₜᵣ)
        tₖᴸ⁻ = boundariesₛₜᵣ[k].tᴸ⁻
        tₖᵁ⁻ = boundariesₛₜᵣ[k].tᵁ⁻
        tₖᴸ⁺ = boundariesₛₜᵣ[k].tᴸ⁺
        tₖᵁ⁺ = boundariesₛₜᵣ[k].tᵁ⁺

        display(L"""
            t_{%$(k)}^{\text{L}-} = %$(tₖᴸ⁻), ~~
            t_{%$(k)}^{\text{U}-} = %$(tₖᵁ⁻), ~~
            t_{%$(k)}^{\text{L}+} = %$(tₖᴸ⁺), ~~
            t_{%$(k)}^{\text{U}+} = %$(tₖᵁ⁺)
        """)
    end
end

ConcatImportanceMethodConcatLPResultLaTeXString = @NamedTuple{
    t⁻::String, t⁺::String,
    W::String,
    vᴸ⁻::String, vᵁ⁻::String,
    vᴸ⁺::String, vᵁ⁺::String,
    εᴸ::String, εᵁ::String,
    }

function concatImportanceMethodConcatLPResultLaTeXString(
        result::ConcatImportanceMethodConcatLPResult{T}
        )::ConcatImportanceMethodConcatLPResultLaTeXString where {T <: Real}
    m = length(result.t⁻)
    n = length(result.wᴸ)

    t⁻ = "\\begin{bmatrix}"; t⁺ = "\\begin{bmatrix}"
    for k = 1:m
        t⁻ *= " $(string(round(result.t⁻[k], digits=3))) "
        t⁺ *= " $(string(round(result.t⁺[k], digits=3))) "
        if k != m
            t⁻ *= " \\\\ "; t⁺ *= " \\\\ "
        end
    end
    t⁻ *= "\\end{bmatrix}"; t⁺ *= "\\end{bmatrix}"

    W = "\\begin{bmatrix}"
    vᴸ⁻ = "\\begin{bmatrix}"; vᵁ⁻ = "\\begin{bmatrix}"
    vᴸ⁺ = "\\begin{bmatrix}"; vᵁ⁺ = "\\begin{bmatrix}"
    εᴸ = "\\begin{bmatrix}"; εᵁ = "\\begin{bmatrix}"
    for i = 1:n
        wᴸ = string(round(result.wᴸ[i], digits=3))
        wᵁ = string(round(result.wᵁ[i], digits=3))
        W *= "\\left[ $(wᴸ), $(wᵁ) \\right]"
        vᴸ⁻ *= " $(string(round(result.vᴸ⁻[i], digits=3))) "
        vᵁ⁻ *= " $(string(round(result.vᵁ⁻[i], digits=3))) "
        vᴸ⁺ *= " $(string(round(result.vᴸ⁺[i], digits=3))) "
        vᵁ⁺ *= " $(string(round(result.vᵁ⁺[i], digits=3))) "
        εᴸ *= " $(string(round(result.εᴸ[i], digits=3))) "
        εᵁ *= " $(string(round(result.εᵁ[i], digits=3))) "
        if i != n
            W *= " \\\\ "
            vᴸ⁻ *= " \\\\ "; vᵁ⁻ *= " \\\\ "
            vᴸ⁺ *= " \\\\ "; vᵁ⁺ *= " \\\\ "
            εᴸ *= " \\\\ "; εᵁ *= " \\\\ "
        end
    end
    W *= "\\end{bmatrix}"
    vᴸ⁻ *= "\\end{bmatrix}"; vᵁ⁻ *= "\\end{bmatrix}"
    vᴸ⁺ *= "\\end{bmatrix}"; vᵁ⁺ *= "\\end{bmatrix}"
    εᴸ *= "\\end{bmatrix}"; εᵁ *= "\\end{bmatrix}"

    return (
        t⁻=t⁻, t⁺=t⁺, W=W,
        vᴸ⁻=vᴸ⁻, vᵁ⁻=vᵁ⁻, vᴸ⁺=vᴸ⁺, vᵁ⁺=vᵁ⁺,
        εᴸ=εᴸ, εᵁ=εᵁ
    )
end

function displayConcatImportanceMethodConcatLPResult(
        result::ConcatImportanceMethodConcatLPResult{T}
        ) where {T <: Real}
    resultₛₜᵣ = concatImportanceMethodConcatLPResultLaTeXString(result)
    
    display(L"""
        t^- = %$(resultₛₜᵣ.t⁻), ~~
        t^+ = %$(resultₛₜᵣ.t⁺)
        """)
    display(L"W = %$(resultₛₜᵣ.W)")
    display(L"""
        v^{\text{L}-} = %$(resultₛₜᵣ.vᴸ⁻), ~~
        v^{\text{U}-} = %$(resultₛₜᵣ.vᵁ⁻), ~~
        v^{\text{L}+} = %$(resultₛₜᵣ.vᴸ⁺), ~~
        v^{\text{U}+} = %$(resultₛₜᵣ.vᵁ⁺)
        """)
    display(L"""
        \varepsilon^\text{L} = %$(resultₛₜᵣ.εᴸ), ~~
        \varepsilon^\text{U} = %$(resultₛₜᵣ.εᵁ)
        """)
end
