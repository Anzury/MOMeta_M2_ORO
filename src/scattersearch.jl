using CSV

include("model.jl")
include("plots.jl")
include("grasp.jl")
include("tabu.jl")
include("pathrelinking.jl")
include("qualitymeasure.jl")

# Run scatter search on all instances in path and save the results in savepath.
# Return a dictionary with the results.
# The parameters are:
#  - path: path to the instances
#  - savepath: path to save the results
#  - popSize: size of the initial population
#  - α: parameter for the GRASP algorithm
#  - p: parameter for the GRASP algorithm
#  - TL: tabu list size
#  - kp: parameter for the tabu search
#  - verboseLevel: 0 for no verbose, 1 for verbose, 2 for verbose (including
#  the verbose messages from the solver used)
#  - exact: true to solve the exact model
#  - unique: true to remove duplicate solutions
#  - timeout: timeout for the exact solver
function ScatterSearch(path::String, savepath::String="", popSize::Int64=10; α::Float64=0.7,
                       p::Float64=0.4, TL::Int64=7, kp::Float64=0.4, verboseLevel::Int64=0,
                       exact=true, unique=false, timeout=typemax(Int64), exactSolSets=nothing)
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

  iter = 1
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
    push!(instime, @elapsed YPN = pathrelinking(ins, imprPop)) # ***
    for sol in YPN
      push!(combiZ1, sol[1])
      push!(combiZ2, sol[2])
    end
    solSet = SolutionSet(ins, length(YPN), YPN)

    if verbose && exact
      println("Solving exact...")
    end
    solSetExact, Q, D = nothing, -1, -1

    if exact
      if isnothing(exactSolSets)
        push!(instime, @elapsed solSetExact = solveExact(path*file, verbose=verbose2, timeout=timeout, unique=unique)) # ***
      else
        texact, solSetExact = exactSolSets[iter]
        push!(instime, texact)
      end
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
    display(plotSolution(solSetExact))
    savefig(savepath*"/map/"*insname*"_exact.png")
    plotAllSolution(solSetExact)
    savefig(savepath*"/map/all_"*insname*"_exact.png")
    display(plotSolution(solSet))
    savefig(savepath*"/map/"*insname*"_SS.png")
    plotAllSolution(solSet)
    savefig(savepath*"/map/all_"*insname*"_SS.png")

    println("----------------------------------\n")
    iter += 1
  end

  return results
end

# Run Scatter Search numIter times and return the average solving time, the average number of solutions
# in YPN, the average quality measure Q and the average distance D for each instance.
# Store the results in a CSV file.
function runScatterSearch(nIter::Int64, path::String, savepath::String="", popSize::Int64=10; α::Float64=0.7,
                       p::Float64=0.6, TL::Int64=7, kp::Float64=0.4, verboseLevel=0, unique=false,
                       timeout=typemax(Int64), exact=true, exactSolSets=nothing)
  if verboseLevel > 0
    println("Scatter Search")
    println("popSize = "*string(popSize))
    println("α = "*string(α))
    println("p = "*string(p))
    println("TL = "*string(TL))
    println("kp = "*string(kp))
    println("verbose level = ", verboseLevel)
    println("unique = "*string(unique))
    println("timeout = "*string(timeout)*"s")
    println("exact = "*string(exact))
    println("\n")
  end
  results = Dict{String, Tuple{Vector{Float64}, Float64, Float64, Float64}}()
  for i=1:nIter
    if verboseLevel > 0
      println("Iteration "*string(i))
    end
    r = ScatterSearch(path, savepath, popSize, α=α, p=p, TL=TL, kp=kp, verboseLevel=-verboseLevel, unique=unique,
                      timeout=timeout, exact=exact, exactSolSets=exactSolSets)

    for (k, v) in r
      if haskey(results, k)
        results[k] = (results[k][1] + v[1], results[k][2] + v[2], results[k][3] + v[3], results[k][4] + v[4])
      else
        results[k] = v
      end
    end
  end

  for (k, v) in results
    results[k] = (v[1]/nIter, v[2]/nIter, v[3]/nIter, v[4]/nIter)
  end

  # sort results in lexicographic order of the instane name
  results = sort(collect(results), by = x -> x[1])
  instances, tGRASP, tTS, tPR, tEX, nSol, Q, D = [], [], [], [], [], [], [], []
  for (k, v) in results
    push!(instances, k)
    push!(tGRASP, round(v[1][1], digits=3))
    push!(tTS, round(v[1][2], digits=3))
    push!(tPR, round(v[1][3], digits=3))
    push!(tEX, round(v[1][4], digits=3))
    push!(nSol, round(v[2], digits=3))
    push!(Q, round(v[3], digits=3))
    push!(D, round(v[4], digits=3))
  end
  results = DataFrame(instance=instances, tGRASP=tGRASP, tTS=tTS, tPR=tPR, tEX=tEX, nSol=nSol, Q=Q, D=D)
  strf = f -> replace(string(round(f, digits=3)), "." => ",")
  S = [string(popSize), strf(α), strf(p), string(TL), strf(kp)]
  fname = savepath*"results_"*S[1]*"_"*S[2]*"_"*S[3]*"_"*S[4]*"_"*S[5]*".csv"
  CSV.write(fname, results)

  return results
end

# Run runScatterSearch with different values for the parameters. For each
# parameters we select uniformly nVal values between the parameters min and max.
# The parameters are:
# - nVal: number of values between min and max
# - nIter: number of iterations for each parameter
# - path: path to the instances
# - savepath: path to save the results
# - popSize: size of the initial population
# - p: parameter for GRASP
# - kp: parameter for tabu search
# - αmin: minimum value for α
# - αmax: maximum value for α
# - TLmin: minimum value for TL
# - TLmax: maximum value for TL
# - verboseLevel: 0 for no verbose, 1 for verbose, 2 for verbose (including
# the verbose messages from the solver used)
# - exact: true to solve the exact model
# - unique: true to remove duplicate solutions
# - timeout: timeout for the exact solver
function runScatterSearchBattery(nVal::Int64, nIter::Int64, path::String, savepath::String, popSize::Int64=10,
                                 p::Float64=0.4, kp::Float64=0.4; αmin::Float64=0.0, αmax::Float64=1.0, TLmin::Int64=3,
                                 TLmax::Int64=14, verboseLevel::Int64=0, exact=true, unique=false, timeout=typemax(Int64))
  if verboseLevel > 0
    println("Scatter Search Tuning")
    println("nVal = "*string(nVal))
    println("nIter = "*string(nIter))
    println("popSize = "*string(popSize))
    println("p = "*string(p))
    println("kp = "*string(kp))
    println("αmin = "*string(αmin))
    println("αmax = "*string(αmax))
    println("TLmin = "*string(TLmin))
    println("TLmax = "*string(TLmax))
    println("verbose level = ", verboseLevel)
    println("unique = "*string(unique))
    println("timeout = "*string(timeout)*"s")
    println("exact = "*string(exact))
    println("Number of configurations to be tested: "*string(nVal^2))
    println("\n")
  end
  αs = range(αmin, αmax, length=nVal)
  TLs = range(TLmin, TLmax, length=nVal)
  exactSolSets = []

  # solve exact and store it in exactSolSets
  results = Dict{String, Vector{Point}}()
  labels = []
  for file in readdir(path)
    if file[1] == '.'
      continue
    end

    t = @elapsed solSetExact = solveExact(path*file, verbose=false, timeout=timeout, unique=unique)
    push!(exactSolSets, (t, solSetExact))
    results[getInstanceName(file)] = []
  end

  for α in αs
    for TL in TLs
      α = round(α, digits=3)
      TL = Int(floor(TL))
      S = [string(popSize), string(α), string(p), string(TL), string(kp)]
      if verboseLevel > 0
        println("Solving with parameters: popSize = "*S[1]*", α = "*S[2]*", p = "*S[3]*", TL = "*S[4]*", kp = "*S[5])
      end
      r = runScatterSearch(nIter, path, savepath, popSize, α=α, p=p, TL=TL, kp=kp, verboseLevel=verboseLevel,
                       unique=unique, timeout=timeout, exact=exact, exactSolSets=exactSolSets)

      for j=1:nrow(r)
        insname = r[j, 1]
        Q = r[j, 7]
        D = r[j, 8]
        push!(results[insname], (Q, D))
      end

      push!(labels, "α = "*S[2]*", TL = "*S[4]*"")
    end
  end

  if !isdir(savepath*"/tuning")
    mkdir(savepath*"/tuning")
  end

  # Plot the results for each instance 
  for (k, L) in results
    iter = 1
    plot(title="Quality measures "*k, xlabel=L"Q_{mean}", ylabel=L"D_{mean}")
    palette(:glasbey_bw_minc_20_hue_150_280_n256)
    for (Q, D) in L
      scatter!([Q], [D], label=labels[iter])
      iter += 1
    end
    savefig(savepath*"/tuning/"*k*".png")
  end
end
