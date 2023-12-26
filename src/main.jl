include("model.jl")
include("plots.jl")
include("grasp.jl")
include("tabu.jl")

function solveAll(path::String, savepath::String=""; verbose=false, timeout=typemax(Int64), unique=false)
  for file in readdir(path)
    if file[1] == '.'
      continue
    end
    insname = getInstanceName(file)
    println("Solving "*insname*"...")
    solSet = solveExact(path*file, verbose=verbose, timeout=timeout, unique=unique)

    if solSet.result_count > 0
      plotYN(solSet)
      savefig(savepath*"/exact/YN/"*insname*".png")
      display(plotSolution(solSet))
      savefig(savepath*"/exact/map/"*insname*".png")
      plotAllSolution(solSet)
      savefig(savepath*"/exact/map/all_"*insname*".png")
    end
  end
end

function testTabu(path::String, savepath::String="", popSize::Int64=10, α::Float64=0.7, TL=7)
  for file in readdir(path)
    if file[1] == '.'
      continue
    end
    insname = getInstanceName(file)
    println("Solving "*insname*"...")

    ins, initPop = createPopulationGRASP(path*file, popSize, α=α)
    initZ1, initZ2 = [], []
    for sol in initPop
      push!(initZ1, sol.z1)
      push!(initZ2, sol.z2)
    end
    plotYN(initZ1, initZ2, lb="GRASP", ms=:+)

    YPN = initSkipList()
    imprZ1, imprZ2 = [], []
    _, solutions = tabuImprovement(ins, initPop, TL=TL, k=0.4*ins.nLvl1)
    for sol in solutions
      push!(YPN, Point(sol.z1, sol.z2))
      push!(imprZ1, sol.z1)
      push!(imprZ2, sol.z2)
    end

    solSet = SolutionSet(ins, length(YPN), YPN, solutions)
    solSetExact = solveExact(path*file, verbose=false)

    plotYN!(imprZ1, imprZ2, lb="Tabu", ms=:x)
    savefig(savepath*"/tabu/obj/"*insname*".png")
    plotYN(solSet, lb=L"Y_{PN}", ms=:x)
    plotYN!(solSetExact)
    savefig(savepath*"/tabu/YN/"*insname*".png")
    display(plotSolution(solSet))
    savefig(savepath*"/tabu/map/"*insname*".png")
    plotAllSolution(solSet)
    savefig(savepath*"/tabu/map/all_"*insname*".png")
  end
end

# solveAll("instances/txt/", "results/txt/")
# solveAll("instances/geojson/", "results/geojson/")

testTabu("instances/txt/", "results/txt/", 100, 0.7, 12)
