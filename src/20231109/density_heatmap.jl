using Plots
pyplot()

function calculate_nbins(range::Tuple{T, T}, bin_size::T) where T <: Real
    return ceil(Int, (range[2] - range[1]) / bin_size)
end

function density_heatmap(
        x::Vector{T},
        y::Vector{T};
        xlabel::AbstractString = "",
        ylabel::AbstractString = "",
        xbinsize::T = (maximum(x) - minimum(x)) / 100,
        ybinsize::T = (maximum(y) - minimum(y)) / 100,
        clims::Tuple{T,T} = (minimum(x), maximum(y)*1.1),
        xlims::Tuple{T,T} = (minimum(x), maximum(x)*1.1),
        ylims::Tuple{T,T} = (minimum(y), maximum(y)*1.1)
        ) where T <: Real
    if length(x) != length(y)
        throw(DimensionMismatch("x and y must have the same length."))
    end

    nbins_x = calculate_nbins(xlims, xbinsize)
    nbins_y = calculate_nbins(ylims, ybinsize)

    return histogram2d(x, y, nbins=(nbins_x, nbins_y),
        xlabel=xlabel,
        ylabel=ylabel,
        color=:viridis,
        colorbar=true,
        clims=clims,
        xlims=xlims,
        ylims=ylims,
    )
end
