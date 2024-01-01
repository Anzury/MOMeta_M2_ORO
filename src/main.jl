include("scattersearch.jl")

# solveAll("instances/txt/", "results/txt/")
# solveAll("instances/geojson/", "results/geojson/")

ScatterSearch("instances/txt/", "results/txt/", 200, α=0.7, verboseLevel=1)
# ScatterSearch("instances/geojson/", "results/geojson/", 100, α=0.7, exact=false, verboseLevel=1)

# runScatterSearch(1, "instances/txt/", "results/txt/", 200, α=0.7, verboseLevel=1)

# runScatterSearchBattery(4, 5, "instances/txt/", "results/txt/", 10, verboseLevel=1)
# runScatterSearchBattery(4, 5, "instances/txt/", "results/txt/", 50, verboseLevel=1)
# runScatterSearchBattery(4, 5, "instances/txt/", "results/txt/", 100, verboseLevel=1)
# runScatterSearchBattery(4, 5, "instances/txt/", "results/txt/", 200, verboseLevel=1)
