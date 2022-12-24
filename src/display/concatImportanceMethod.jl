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

    # W⁻ = ∅
    if any(isinf.(result.wₖᴸ⁻)) || any(isinf.(result.wₖᵁ⁻))
        wₖᴸ⁺ = "\\begin{bmatrix}"; wₖᵁ⁺ = "\\begin{bmatrix}"
        for i = 1:n
            wₖᴸ⁺ *= " $(string(round(result.wₖᴸ⁺[i], digits=3))) "
            wₖᵁ⁺ *= " $(string(round(result.wₖᵁ⁺[i], digits=3))) "
            if i != n
                wₖᴸ⁺ *= " \\\\ "; wₖᵁ⁺ *= " \\\\ "
            end
        end
        wₖᴸ⁺ *= "\\end{bmatrix}"; wₖᵁ⁺ *= "\\end{bmatrix}"
    
        return (wₖᴸ⁻="", wₖᵁ⁻="", wₖᴸ⁺=wₖᴸ⁺, wₖᵁ⁺=wₖᵁ⁺)
    end

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

        if resultₛₜᵣ.wₖᴸ⁻ == "" || resultₛₜᵣ.wₖᵁ⁻ == ""
            display(L"W_%$(k)^- = \emptyset")
            display(L"""
                w_{%$(k)}^{\text{L}+} = %$(resultₛₜᵣ.wₖᴸ⁺), ~~
                w_{%$(k)}^{\text{U}+} = %$(resultₛₜᵣ.wₖᵁ⁺)
            """)
        else
            display(L"""
                w_{%$(k)}^{\text{L}-} = %$(resultₛₜᵣ.wₖᴸ⁻), ~~
                w_{%$(k)}^{\text{U}-} = %$(resultₛₜᵣ.wₖᵁ⁻), ~~
                w_{%$(k)}^{\text{L}+} = %$(resultₛₜᵣ.wₖᴸ⁺), ~~
                w_{%$(k)}^{\text{U}+} = %$(resultₛₜᵣ.wₖᵁ⁺)
            """)
        end
    end
end

ConcatImportanceMethodTBoundariesLaTeXString = @NamedTuple{
    tₖᴸ⁻::String, tₖᵁ⁻::String, tₖᴸ⁺::String, tₖᵁ⁺::String
    }

function concatImportanceMethodTBoundariesLaTeXString(
        boundaries::ConcatImportanceMethodTBoundaries{T}
        )::ConcatImportanceMethodTBoundariesLaTeXString where {T <: Real}
    if isinf(boundaries.tₖᴸ⁻) || isinf(boundaries.tₖᵁ⁻)
        tₖᴸ⁺ = string(round(boundaries.tₖᴸ⁺, digits=3))
        tₖᵁ⁺ = string(round(boundaries.tₖᵁ⁺, digits=3))
        return (tₖᴸ⁻="", tₖᵁ⁻="", tₖᴸ⁺=tₖᴸ⁺, tₖᵁ⁺=tₖᵁ⁺)
    else
        tₖᴸ⁻ = string(round(boundaries.tₖᴸ⁻, digits=3))
        tₖᵁ⁻ = string(round(boundaries.tₖᵁ⁻, digits=3))
        tₖᴸ⁺ = string(round(boundaries.tₖᴸ⁺, digits=3))
        tₖᵁ⁺ = string(round(boundaries.tₖᵁ⁺, digits=3))
        return (tₖᴸ⁻=tₖᴸ⁻, tₖᵁ⁻=tₖᵁ⁻, tₖᴸ⁺=tₖᴸ⁺, tₖᵁ⁺=tₖᵁ⁺)
    end
end

function displayConcatImportanceMethodTBoundaries(
        boundaries::ConcatImportanceMethodTBoundaries{T},
        k::Integer) where {T <: Real}
    boundariesₛₜᵣ = concatImportanceMethodTBoundariesLaTeXString(boundaries)

    tₖᴸ⁻ = boundariesₛₜᵣ.tₖᴸ⁻
    tₖᵁ⁻ = boundariesₛₜᵣ.tₖᵁ⁻
    tₖᴸ⁺ = boundariesₛₜᵣ.tₖᴸ⁺
    tₖᵁ⁺ = boundariesₛₜᵣ.tₖᵁ⁺

    if tₖᴸ⁻  == "" || tₖᵁ⁻ == ""
        display(L"""
            t_{%$(k)}^{\text{L}+} = %$(tₖᴸ⁺), ~~
            t_{%$(k)}^{\text{U}+} = %$(tₖᵁ⁺)
        """)
    else
        display(L"""
            t_{%$(k)}^{\text{L}-} = %$(tₖᴸ⁻), ~~
            t_{%$(k)}^{\text{U}-} = %$(tₖᵁ⁻), ~~
            t_{%$(k)}^{\text{L}+} = %$(tₖᴸ⁺), ~~
            t_{%$(k)}^{\text{U}+} = %$(tₖᵁ⁺)
        """)
    end
end

function displayConcatImportanceMethodTBoundaries(
        boundaries::AbstractArray{ConcatImportanceMethodTBoundaries{T}}) where {T <: Real}
    for k = eachindex(boundaries)
        displayConcatImportanceMethodTBoundaries(boundaries[k], k)
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
