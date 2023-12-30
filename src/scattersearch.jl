include("model.jl")
include("plots.jl")
include("grasp.jl")
include("tabu.jl")
include("pathrelinking.jl")
include("qualitymeasure.jl")

function ScatterSearch(path::String, savepath::String="", popSize::Int64=10; α::Float64=0.7,
                       p::Float64=0.4, TL::Int64=7, kp::Float64=0.4, verboseLevel::Int64=0,
                       exact=true, unique=false, timeout=typemax(Int64))
  verbose, verbose2 = verboseLevel^2 > 0, verboseLevel == 2 || verboseLevel == -2
  if verbose && verboseLevel >= 0
    println("Scatter Search")
    println("popSize = "*string(popSize))
    println("α = "*string(α))
    println("p = "*string(p))
    println("TL = "*string(TL))
    println("kp = "*string(kp))
    println("verbose = "*string(verbose))
    println("unique = "*string(unique))
    println("timeout = "*string(timeout)*"s")
    println("exact = "*string(exact))
    println("\n")
  end

  if !isdir(savepath*"/YN")
    mkdir(savepath*"/YN")
  end
  if !isdir(savepath*"/map")
    mkdir(savepath*"/map")
  end
  if !isdir(savepath*"/obj")
    mkdir(savepath*"/obj")
  end

  results = Dict{String, Tuple{Vector{Float64}, Float64, Float64, Float64}}()
  for file in readdir(path)
    if file[1] == '.'
      continue
    end
    instime::Vector{Float64} = []
    insname = getInstanceName(file)
    println("----------------------------------")
    println("Solving "*insname*"\n")

    if verbose
      println("Creating initial population...")
    end
    push!(instime, @elapsed ins, initPop = createPopulationGRASP(path*file, popSize, α=α, p=p)) # ***
    initZ1, initZ2 = [], []
    for sol in initPop
      push!(initZ1, sol.z1)
      push!(initZ2, sol.z2)
    end
    plotYN(initZ1, initZ2, lb="GRASP", ms=:+)

    if verbose
      println("Improving initial population...")
    end
    imprZ1, imprZ2 = [], []
    push!(instime, @elapsed _, imprPop = tabuImprovement(ins, initPop, TL=TL, k=kp*ins.nLvl1)) # ***
    for sol in imprPop
      push!(imprZ1, sol.z1)
      push!(imprZ2, sol.z2)
    end

    if verbose
      println("Path relinking...")
    end
    combiZ1, combiZ2 = [], []
    push!(instime, @elapsed YPN, combiPop = pathrelinking(ins, imprPop)) # ***
    for sol in combiPop
      push!(combiZ1, sol.z1)
      push!(combiZ2, sol.z2)
    end
    finalPop = [imprPop; combiPop]
    # remove all solutions whose objective values are not in YPN
    finalPop = filter(sol -> (sol.z1, sol.z2) in YPN, finalPop)
    # order finalPop as in YPN
    finalPop = sort!(finalPop, by = x -> (x.z1, x.z2))
    solSet = SolutionSet(ins, length(YPN), YPN, finalPop)

    if verbose && exact
      println("Solving exact...")
    end
    solSetExact, Q, D = nothing, -1, -1
    if exact
      push!(instime, @elapsed solSetExact = solveExact(path*file, verbose=verbose2, timeout=timeout, unique=unique)) # ***
      areaYPN = computeArea(solSetExact.YN, YPN)
      areaexact = computeAreaExact(solSetExact.YN, YPN)
      Q = areaexact != 0 ? areaYPN/areaexact : -0
      Q = isnan(Q) ? -0 : Q*100
      D = computeIdealsDistance(solSetExact.YN, YPN)
    else
      push!(instime, -1.0)
    end
    results[insname] = (instime, length(YPN)+0.0, Q, D)

    plotYN!(imprZ1, imprZ2, lb="Tabu", ms=:x)
    plotYN!(combiZ1, combiZ2, lb="Path relinking", ms=:star4)
    if exact
      plotYN!(solSetExact)
    end
    savefig(savepath*"/obj/"*insname*".png")
    plotYN(solSet, lb=L"Y_{PN}", ms=:x)
    if exact
      plotYN!(solSetExact)
    end
    savefig(savepath*"/YN/"*insname*".png")
    display(plotSolution(solSet))
    savefig(savepath*"/map/"*insname*".png")
    plotAllSolution(solSet)
    savefig(savepath*"/map/all_"*insname*".png")

    println("----------------------------------\n")
  end

  return results
end

