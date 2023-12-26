include("parser.jl")

# Get greedy function value for concentrator j w.r.t 1st obj function
function g1(ins::Instance, j::Int64)
  # g1 is computed as the sum of ins.C[i,j] for all i terminals and lvl1
  # concentrator j + the sum of ins.B[j,k] for all k level 2 concentrators
  # since we minimize the first objective function we return 1/g1
  return 1/sum(ins.C[i,j] for i in 1:ins.nTerm) + sum(ins.B[j,k] for k in 1:ins.nLvl2)
end

# Get greedy function value for concentrator j w.r.t 2nd obj function
function g2(ins::Instance, j::Int64)
  # g2 is computed as the maximum of ins.C[i,j] for all i terminals and lvl1
  # concentrator j
  # since we minimize the second objective function we return 1/g2
  return 1/maximum([ins.C[i,j] for i in 1:ins.nTerm])
end

# GRASP method : returns a greedy solution w.r.t to a greedy function, uses
# a parameter α and p representing the proportion of concentrators to
# locate
function GRASP(ins::Instance, g::Function; α::Float64=0.7, p::Float64=0.4)
  x = zeros(Int64, ins.nTerm, ins.nLvl1)
  y = zeros(Int64, ins.nLvl1, ins.nLvl2)
  z = zeros(Int64, ins.nLvl2)
  U = [g(ins, j) for j in 1:ins.nLvl1]
  U_index = sortperm(U)
  # the number of concentrator to locate
  p = Int(ceil(p*ins.nLvl1))

  j = rand(1:ins.nLvl1) # we first select a lvl1 concentrator randomly 
  lvl1 = [j] # store potential lvl1 concentrators

  # we remove j from the candidate list by setting U[j] = -1
  # From here on, to check if a concentrator is a candidate we check if U[j] != -1
  U[j] = -1

  # let gmax and gmin be the max and min of U_index (to find the min we check
  # the beginning of U_index and to find the max we check the end of U_index,
  # the max/min is the first element in U_index such that U[U_index[i]] != -1)
  function getgvalues(localU, localU_index)
    localgmin, localgmax = -1, -1
    for i in 1:ins.nLvl1
      if U[U_index[i]] != -1
        localgmin = U[U_index[i]]
        break
      end
    end

    for i in ins.nLvl1:-1:1
      if U[U_index[i]] != -1
        localgmax = U[U_index[i]]
        break
      end
    end

    return localgmin, localgmax
  end

  gmin, gmax = getgvalues(U, U_index)

  # we compute L the limit such as L = gmax - α * (gmax - gmin)
  L = gmax - α * (gmax - gmin)

  # Select lvl1 concentrator iteratively until p concentrators are selected
  for i=1:p
    # construct RCL
    RCL = [j for j in 1:ins.nLvl1 if U[j] >= L]
    if RCL == []
      break
    end
    j = rand(RCL) # we select a lvl1 concentrator randomly in RCL
    push!(lvl1, j) # we add j to the lvl1 concentrator list
    U[j] = -1 # we remove j from the candidate list
    gmin, gmax = getgvalues(U, U_index) # we update gmin and gmax
    L = gmax - α * (gmax - gmin) # we update L
  end

  # Connect each terminal to its two closest lvl1 concentrators
  for i in 1:ins.nTerm
    # we select the two closest lvl1 concentrators to terminal i
    D = [ins.C[i,j] for j in lvl1]
    j1, j2, d1, d2 = -1, -1, typemax(Int64), typemax(Int64)
    for k=1:length(D)
      if D[k] < d1
        j1, d1 = k, D[k]
      end
    end
    for k=1:length(D)
      if D[k] < d2 && k != j1
        j2, d2 = lvl1[k], D[k]
      end
    end
    j1 = lvl1[j1]

    x[i,j1] = 1
    x[i,j2] = 1
  end

  # Now that the lvl1 concentrators are connected to terminals, we connect them
  # to lvl2 concentrators such that each lvl1 concentrator is connected to
  # exactly one lvl2 concentrator. If a lvl1 concentrator is not connected to
  # any terminal, we do not connect it to any lvl2 concentrator.
  for j in 1:ins.nLvl1
    if sum(x[i,j] for i in 1:ins.nTerm) > 0
      # we select a lvl2 concentrator k such as ins.B[j,k] is minimum
      k = sortperm([ins.B[j,k] for k in 1:ins.nLvl2])[1]
      y[j,k] = 1
      z[k] = 1
    end
  end

  return x, y, z
end

# Create a population of size popSize using GRASP
function createPopulationGRASP(path::String, popSize::Int64; α::Float64=0.7)
  # Load instance
  ins::Instance = loadInstance(path)

  x = zeros(Int64, ins.nTerm, ins.nLvl1)
  y = zeros(Int64, ins.nLvl1, ins.nLvl2)
  z = zeros(Int64, ins.nLvl2)
  population::Vector{Solution} = []

  # Create population
  halfpop = Int(ceil(popSize/2))
  for _ in 1:halfpop # create popSize/2 solutions with g1
    x, y, z = GRASP(ins, g1, α=α)
    z1, z2 = getZ1(ins, x, y, z), getZ2(ins, x)
    push!(population, Solution(x, y, z, z1, z2))
  end

  for _ in halfpop+1:popSize # create popSize/2 solutions with g2
    x, y, z = GRASP(ins, g2, α=α)
    z1, z2 = getZ1(ins, x, y, z), getZ2(ins, x)
    push!(population, Solution(x, y, z, z1, z2))
  end

  return ins, population
end
