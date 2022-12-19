using IntervalArithmetic

include("../intervalPCM/intervalPCM.jl")

"""
区間ベクトルを LaTeX 形式にする  
`L"intervalVectorLaTeXString(W)"` とすると LaTeX 形式で表示できる
"""
function intervalVectorLaTeXString(
        W::Vector{Interval{T}}
        )::String where {T <: Real}
    n = length(W)

    # 各成分の LaTeX 表記を入れる
    Wₛₜᵣ = fill("", n)
    for i = 1:n
        # (i,j)成分の両端
        wᵢᴸ = string(round(W[i].lo, digits=3))
        wᵢᵁ = string(round(W[i].hi, digits=3))
        Wₛₜᵣ[i] = "\\left[ $(wᵢᴸ), $(wᵢᵁ) \\right]"
    end

    str = "\\begin{bmatrix} "
    for i = 1:n
        str = str * "$(Wₛₜᵣ[i])"
        if i != n
            str = str * " \\\\ "
        end
    end
    str = str * " \\end{bmatrix}"

    return str
end

"""
区間行列を LaTeX 形式にする  
`L"intervalMatrixLaTeXString(A)"` とすると LaTeX 形式で表示できる
"""
function intervalMatrixLaTeXString(
        A::Matrix{Interval{T}}
        )::String where {T <: Real}
    m, n = size(A)

    # 各成分の LaTeX 表記を入れる
    Aₛₜᵣ = fill("", m, n)
    for i = 1:m, j = 1:n
        # (i,j)成分の両端
        aᵢⱼᴸ = string(round(A[i,j].lo, digits=3))
        aᵢⱼᵁ = string(round(A[i,j].hi, digits=3))
        Aₛₜᵣ[i,j] = "\\left[ $(aᵢⱼᴸ), $(aᵢⱼᵁ) \\right]"
    end

    str = "\\begin{bmatrix} "
    for i = 1:m, j = 1:n
        str = str * "$(Aₛₜᵣ[i,j])"
        if j == n
            if i != m
                str = str * " \\\\ "
            end
        else
            str = str * " & "
        end
    end
    str = str * " \\end{bmatrix}"

    return str
end
