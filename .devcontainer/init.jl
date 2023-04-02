using Pkg

function main()
    Pkg.add([
        "JuMP",
        "HiGHS",
        "IntervalArithmetic",
        "Printf",
        "Latexify",
        "LaTeXStrings",
        "Plots",
        "PyPlot",
        "Distributions"
        ])
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
