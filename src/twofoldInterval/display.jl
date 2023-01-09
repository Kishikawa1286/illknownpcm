include("./twofoldInterval.jl")

"""
二重区間行列を LaTeX 形式にする  
`L"twofoldIntervalMatrixLaTeXString(A)"` とすると LaTeX 形式で表示できる
"""
function twofoldIntervalMatrixLaTeXString(
        A::Matrix{TwofoldInterval{T}}
        )::String where {T <: Real}
    m, n = size(A)

    # 各成分の LaTeX 表記を入れる
    Aₛₜᵣ = fill("", m, n)
    for i = 1:m, j = 1:n
        # (i,j)成分の外の両端
        aᵢⱼᴸ⁺ = lpad(string(round(A[i,j][2].lo, digits=3)), 5)
        aᵢⱼᵁ⁺ = lpad(string(round(A[i,j][2].hi, digits=3)), 5)
        if iscommon(A[i,j][1])
            # (i,j)成分の内の両端
            aᵢⱼᴸ⁻ = lpad(string(round(A[i,j][1].lo, digits=3)), 5)
            aᵢⱼᵁ⁻ = lpad(string(round(A[i,j][1].hi, digits=3)), 5)
            Aₛₜᵣ[i,j] = "\\left[ $(aᵢⱼᴸ⁺), \\left[ $(aᵢⱼᴸ⁻), $(aᵢⱼᵁ⁻) \\right], $(aᵢⱼᵁ⁺) \\right]"
        else
            Aₛₜᵣ[i,j] = "\\left[ $(aᵢⱼᴸ⁺), \\emptyset, $(aᵢⱼᵁ⁺) \\right]"
        end
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
