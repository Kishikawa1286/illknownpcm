using Pkg

function main()
    Pkg.add(
        Pkg.PackageSpec[
            Pkg.PackageSpec(name="Distributions", version="0.25.100"),
            Pkg.PackageSpec(name="HiGHS", version="1.6.0"),
            Pkg.PackageSpec(name="IntervalArithmetic", version="0.21.0"),
            Pkg.PackageSpec(name="JuMP", version="1.14.1"),
            Pkg.PackageSpec(name="LaTeXStrings"),
            Pkg.PackageSpec(name="Latexify"),
            Pkg.PackageSpec(name="Plots", version="1.15.2"),
            Pkg.PackageSpec(name="Printf"),
            Pkg.PackageSpec(name="ProgressMeter"),
            Pkg.PackageSpec(name="PyPlot", version="2.11.2"),
            Pkg.PackageSpec(name="StatsPlots", version="0.14.33")
        ]
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
