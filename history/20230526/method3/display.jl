using LaTeXStrings

include("./method3.jl")

ApproximationLPResultLaTeXString_m3 = @NamedTuple{
    Wₖ⁻::String, Wₖ⁺::String
    }

function approximationLPResultLaTeXString_m3(
        result::ApproximationLPResult_m3{T}
        )::ApproximationLPResultLaTeXString_m3 where {T <: Real}
    n = length(result.wₖᴸ⁻)

    # W⁻ = ∅
    if any(isinf.(result.wₖᴸ⁻)) || any(isinf.(result.wₖᵁ⁻))
        Wₖ⁺ = "\\begin{bmatrix}"
        for i = 1:n
            wₖᴸ⁺ = string(round(result.wₖᴸ⁺[i], digits=3))
            wₖᵁ⁺ = string(round(result.wₖᵁ⁺[i], digits=3))
            Wₖ⁺ *= " \\left[ $(wₖᴸ⁺), $(wₖᵁ⁺) \\right] "
            if i != n
                Wₖ⁺ *= " \\\\ "
            end
        end
        Wₖ⁺ *= "\\end{bmatrix}"
    
        return (Wₖ⁻="", Wₖ⁺=Wₖ⁺)
    end

    Wₖ⁻ = "\\begin{bmatrix}"; Wₖ⁺ = "\\begin{bmatrix}"
    for i = 1:n
        wₖᴸ⁻ = string(round(result.wₖᴸ⁻[i], digits=3))
        wₖᵁ⁻ = string(round(result.wₖᵁ⁻[i], digits=3))
        Wₖ⁻ *= " \\left[ $(wₖᴸ⁻), $(wₖᵁ⁻) \\right] "
        wₖᴸ⁺ = string(round(result.wₖᴸ⁺[i], digits=3))
        wₖᵁ⁺ = string(round(result.wₖᵁ⁺[i], digits=3))
        Wₖ⁺ *= " \\left[ $(wₖᴸ⁺), $(wₖᵁ⁺) \\right] "
        if i != n
            Wₖ⁻ *= " \\\\ "; Wₖ⁺ *= " \\\\ "
        end
    end
    Wₖ⁻ *= "\\end{bmatrix}"; Wₖ⁺ *= "\\end{bmatrix}"

    return (Wₖ⁻=Wₖ⁻, Wₖ⁺=Wₖ⁺)
end

function displayApproximationLPResults_m3(
        results::AbstractArray{ApproximationLPResult_m3{T}}
        ) where {T <: Real}
    for k = eachindex(results)
        resultₛₜᵣ = approximationLPResultLaTeXString_m3(results[k])

        if resultₛₜᵣ.Wₖ⁻ == ""
            display(L"W_%$(k)^- = \\emptyset")
            display(L"W_%$(k)^+ = %$(resultₛₜᵣ.Wₖ⁺)")
        else
            display(L"W_%$(k)^- = %$(resultₛₜᵣ.Wₖ⁻)")
            display(L"W_%$(k)^+ = %$(resultₛₜᵣ.Wₖ⁺)")
        end
    end
end

TBoundariesLaTeXString_m3 = @NamedTuple{
    tₖᴸ⁻::String, tₖᵁ⁻::String, tₖᴸ⁺::String, tₖᵁ⁺::String
    }

function tBoundariesLaTeXString_m3(
        boundaries::TBoundaries_m3{T}
        )::TBoundariesLaTeXString_m3 where {T <: Real}
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

function displayTBoundaries_m3(
        boundaries::TBoundaries_m3{T},
        k::Integer) where {T <: Real}
    boundariesₛₜᵣ = tBoundariesLaTeXString_m3(boundaries)

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

function displayTBoundaries_m3(
        boundaries::AbstractArray{TBoundaries_m3{T}}) where {T <: Real}
    for k = eachindex(boundaries)
        displayTBoundaries_m3(boundaries[k], k)
    end
end

ConcatLPResultLaTeXString_m3 = @NamedTuple{
    t⁻::String, t⁺::String,
    W::String,
    vᴸ⁻::String, vᵁ⁻::String,
    vᴸ⁺::String, vᵁ⁺::String,
    εᴸ::String, εᵁ::String,
    }

function concatLPResultLaTeXString_m3(
        result::ConcatLPResult_m3{T}
        )::ConcatLPResultLaTeXString_m3 where {T <: Real}
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

function displayConcatLPResult_m3(
        result::ConcatLPResult_m3{T}
        ) where {T <: Real}
    resultₛₜᵣ = concatLPResultLaTeXString_m3(result)
    
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
