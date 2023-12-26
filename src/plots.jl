using Plots, LaTeXStrings

function getZ1Z2fromSolSet(s::SolutionSet)
  Z1, Z2 = [], []
  if s.result_count < 1
    return Z1, Z2
  end
  for p in s.YN
    push!(Z1, p[1])
    push!(Z2, p[2])
  end
  return Z1, Z2
end

# function to plot YN
function plotYN(s::SolutionSet; lb=L"Y_N", ms=:dot, st=:scatter)
  Z1, Z2 = getZ1Z2fromSolSet(s)
  plotYN(Z1, Z2, lb=lb, ms=ms, st=st)
end

function plotYN!(s::SolutionSet; lb=L"Y_N", ms=:dot, st=:scatter)
  Z1, Z2 = getZ1Z2fromSolSet(s)
  plotYN!(Z1, Z2, lb=lb, ms=ms, st=st)
end

function plotYN(Z1, Z2; lb=L"Y_N", ms=:dot, st=:scatter)
  plot(Z1, Z2, seriestype=st, marker=(ms, 4), label=lb, xlabel=L"z_1", ylabel=L"z_2", title="Objective space")
end

function plotYN!(Z1, Z2; lb=L"Y_N", ms=:dot, st=:scatter)
  plot!(Z1, Z2, seriestype=st, marker=(ms, 4), label=lb, xlabel=L"z_1", ylabel=L"z_2", title="Objective space")
end

# plot solution obtained after solving the model (by
# default, plots the first nondominated solution)
function plotSolution(s::SolutionSet, r=1)
  if s.result_count < 1
    return
  end
  I = s.instance.nTerm
  J = s.instance.nLvl1
  K = s.instance.nLvl2

  z = (Int(s.YN[r][1]), Int(s.YN[r][2]))
  p = plot(title=string(z), legend=false, showaxis=false)
  
  lvl1term = [[], []]
  for i=1:I
    for j=1:J
      if s.solutions[r].x[i,j] == 1.0
        for t=1:2
          push!(lvl1term[t], [s.instance.terminals[i,t], s.instance.locationLvl1[j,t]])
        end
      end
    end
  end
  plot!(lvl1term[1], lvl1term[2], color=:green, label="")

  lvl1lvl2 = [[], []]
  for j=1:J
    for k=1:K
      if s.solutions[r].y[j,k] == 1.0
        for t=1:2
          push!(lvl1lvl2[t], [s.instance.locationLvl2[k,t], s.instance.locationLvl1[j,t]])
        end
      end
    end
  end
  plot!(lvl1lvl2[1], lvl1lvl2[2], color=:blue, label="")

  scatter!(s.instance.terminals[:, 1], s.instance.terminals[:, 2], color=:red, label="Terminals")
  scatter!(s.instance.locationLvl1[:, 1], s.instance.locationLvl1[:, 2], color=:blue, shape=:rect, label="L1 concentrators")
  scatter!(s.instance.locationLvl2[:, 1], s.instance.locationLvl2[:, 2], color=:black, shape=:rect, label="L2 concentrators")

  return p
end

# plot all nondominated solutions obtained after
# solving the model (or plot only `max` solutions)
function plotAllSolution(solution::SolutionSet; max=30)
  if solution.result_count < 1
    return
  end
  plots = []

  for i=1:max
    if i <= solution.result_count
      push!(plots, plotSolution(solution, i))
    end
  end

  insname = getInstanceName(solution.instance.fname)
  return plot(plots..., size=(1920, 1080), plot_title="All solutions "*insname)
end
