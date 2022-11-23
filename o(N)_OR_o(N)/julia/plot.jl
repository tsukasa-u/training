# import Pkg
# Pkg.add("Plots")

# Pkg.add("PGFPlotsX")

# Pkg.add("UnicodePlots")

# Pkg.add("PyPlot")

# Pkg.add("InspectDR")
# Pkg.add("Gaston")

# Pkg.add("StatsPlots")
# Pkg.add("GraphRecipes")

# using Plots

# gr(size = (300, 300), legend = false)  # provide optional defaults
# pgfplotsx()
# plotly(ticks=:native)                  # plotlyjs for richer saving options
# pyplot()                               # backends are selected with lowercase names
# unicodeplots()

using Plots
x = 1:10; y = rand(10); # These are the plotting data
plot(x, y)










