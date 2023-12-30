using JuMP
import HiGHS, Gurobi
import MultiObjectiveAlgorithms as MOA

include("parser.jl")
include("plots.jl")

function getModelTSUFLP(ins::Instance; solver=Gurobi.Optimizer, timeout=typemax(Int64), unique=false)
  # Model definition (w/ algorithm used)
  model = JuMP.Model(() -> MOA.Optimizer(solver))
  set_attribute(model, MOA.Algorithm(), MOA.EpsilonConstraint())
  set_attribute(model, MOI.TimeLimitSec(), timeout)
  if unique
    set_attribute(model, MOA.SolutionLimit(), 1)
  end

  I = ins.nTerm
  J = ins.nLvl1
  K = ins.nLvl2

  # Variables definition
  @variable(model, x[1:I,1:J], Bin)
  @variable(model, y[1:J,1:K], Bin)
  @variable(model, z[1:K], Bin)
  @variable(model, M >= 0)

  # Expressions definition
  @expression(model, z1_expr, sum(ins.C[i,j]*x[i,j] for j=1:J for i=1:I)
                              + sum(ins.B[j,k]*y[j,k] for k=1:K for j=1:J)
                              + sum(ins.S[k]*z[k] for k=1:K)) 
  @expression(model, z2_expr, M) # pull objective

  # Objectives definitions
  @objective(model, Min, [z1_expr, z2_expr])

  # Constraints definition
  # Each terminal is connected to exactly two concentrators (multiple homing)
  @constraint(model, [i=1:I], sum(x[i,j] for j=1:J) == 2)
  # lvl1 and lvl2 concentrators connection
  @constraint(model, [i=1:I, j=1:J], x[i,j] <= sum(y[j,k] for k=1:K))
  # Each lvl1 concentrators are connected to at most one lvl2 concentrator
  @constraint(model, [j=1:J, k=1:K], y[j,k] <= z[k])
  @constraint(model, [j=1:J], sum(y[j,k] for k=1:K) <= 1)
  # Multiple homing constraint to have M >= max(C[i,j]*x[i,j])
  @constraint(model, [i=1:I, j=1:J], M >= ins.C[i,j]*x[i,j])

  return model
end

function solveExact(path::String; solver=Gurobi.Optimizer, verbose=true, timeout=typemax(Int64), unique=false)
  # Load instance
  instance::Instance = loadInstance(path)

  # Get model
  model = getModelTSUFLP(instance, solver=solver, timeout=timeout, unique=unique)

  # Solve model
  if !verbose
    set_attribute(model, "output_flag", false)
  end
  time
  optimize!(model)

  R = result_count(model)
  solutions::Vector{Solution} = []
  YN = initSkipList()
  x = zeros(Int, instance.nTerm, instance.nLvl1)
  y = zeros(Int, instance.nLvl1, instance.nLvl2)
  z = zeros(Int, instance.nLvl2)
  for r=1:R
    z1, z2 = objective_value(model, result=r)
    push!(YN, Point(z1, z2))

    for i=1:instance.nTerm
      for j=1:instance.nLvl1
        x[i,j] = value(model[:x][i,j], result=r)
        if i == 1
          for k=1:instance.nLvl2
            y[j,k] = value(model[:y][j,k], result=r)
            z[k] = value(model[:z][k], result=r)
          end
        end
      end
    end

    push!(solutions, Solution(x, y, z, z1, z2))
  end

  if verbose
    solution_summary(model)
  end

  return SolutionSet(instance, R, YN, solutions)
end

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
      savefig(savepath*"/YN_"*insname*".png")
      display(plotSolution(solSet))
      savefig(savepath*"/map_"*insname*".png")
      plotAllSolution(solSet)
      savefig(savepath*"/allmap_"*insname*".png")
    end
  end
end

