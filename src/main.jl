using CSV

include("scattersearch.jl")

# Run Scatter Search numIter times and return the average solving time and the average number of solutions
# in YPN. Store the results in a CSV file.
function runScatterSearch(numIter::Int64, path::String, savepath::String="", popSize::Int64=10; α::Float64=0.7,
                       p::Float64=0.4, TL::Int64=7, kp::Float64=0.4, verboseLevel=0, unique=false,
                       timeout=typemax(Int64), exact=true)
  if verboseLevel > 0
    println("Scatter Search Battery test")
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
  for i=1:numIter
    if verboseLevel > 0
      println("Iteration "*string(i))
    end
    r = ScatterSearch(path, savepath, popSize, α=α, p=p, TL=TL, kp=kp, verboseLevel=-verboseLevel, unique=unique,
                      timeout=timeout, exact=exact)

    for (k, v) in r
      if haskey(results, k)
        results[k] = (results[k][1] + v[1], results[k][2] + v[2], results[k][3] + v[3], results[k][4] + v[4])
      else
        results[k] = v
      end
    end
  end

  for (k, v) in results
    results[k] = (v[1]/numIter, v[2]/numIter, v[3]/numIter, v[4]/numIter)
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
  CSV.write(savepath*"results.csv", results)

  return results
end

# solveAll("instances/txt/", "results/txt/")
# solveAll("instances/geojson/", "results/geojson/")

# ScatterSearch("instances/txt/", "results/txt/", 200, α=0.7, verboseLevel=1)
# ScatterSearch("instances/geojson/", "results/geojson/", 100, α=0.7, exact=false, verboseLevel=1)

runScatterSearch(10, "instances/txt/", "results/txt/", 200, α=0.7, verboseLevel=1)
