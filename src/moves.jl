include("neighborhoods.jl")

function getNonTabuSwapMoves(ins, sol, TM, bestValue; firstobj=true)
  id, iter, maxIter = 0, 0, 100
  lvl1 = [j for j in 1:ins.nLvl1 if isconnected(sol.x[:,j])]

  while iter < maxIter && length(lvl1) > 0
    j1, j2 = rand(lvl1), rand(lvl1)
    if j1 != j2 && !((id, j1, j2) in TM) && !((id, j2, j1) in TM)
      return (id, j1, j2)
    elseif j1 != j2
      new_sol = swap(ins, sol, j1, j2)
      val = firstobj ? new_sol.z1 : new_sol.z2
      if val < bestValue
        return (id, j1, j2)
      end
    end
    iter += 1
  end
  return (id, 1, 1)
end

function getNonTabuShiftMoves(ins, sol, TM, bestValue; firstobj=true)
  id, iter, maxIter = 1, 0, 100
  lvl1 = [j for j in 1:ins.nLvl1 if isconnected(sol.x[:,j])]

  while iter < maxIter && length(lvl1) > 0
    i = -1
    j1 = rand(lvl1)
    j2 = rand(lvl1)
    # select randomly one terminal assigned to j1
    # but not assigned to j2
    for ii in shuffle(1:ins.nTerm)
      if sol.x[ii, j1] == 1 && sol.x[ii, j2] == 0
        i = ii
        break
      end
    end

    if i == -1
      continue
    end

    if j1 != j2 && !((id, i, j1, j2) in TM) && !((id, i, j2, j1) in TM)
      return (id, i, j1, j2)
    elseif j1 != j2
      new_sol = shift(ins, sol, i, j1, j2)
      val = firstobj ? new_sol.z1 : new_sol.z2
      if val < bestValue
        return (id, i, j1, j2)
      end
    end
    iter += 1
  end
  return (id, 1, 1, 1)
end

function getNonTabuDropMoves(ins, sol, TM, bestValue; firstobj=true)
  id, iter, maxIter = 2, 0, 100
  lvl1 = [j for j in 1:ins.nLvl1 if isconnected(sol.x[:,j])]

  while iter < maxIter && length(lvl1) > 0
    j = rand(lvl1)
    if !((id, j) in TM)
      return (id, j)
    else
      new_sol = drop(ins, sol, j)
      val = firstobj ? new_sol.z1 : new_sol.z2
      if val < bestValue
        return (id, j)
      end
    end
    iter += 1
  end
  return (id, 1)
end

function getNonTabuOpenLvl1Moves(ins, sol, TM, bestValue; firstobj=true)
  id, iter, maxIter = 3, 0, 100
  lvl1 = [j for j in 1:ins.nLvl1 if !isconnected(sol.x[:,j])]

  while iter < maxIter && length(lvl1) > 0
    j = rand(lvl1)
    if !((id, j) in TM)
      return (id, j)
    else
      new_sol = openLvl1(ins, sol, j)
      val = firstobj ? new_sol.z1 : new_sol.z2
      if val < bestValue
        return (id, j)
      end
    end
    iter += 1
  end
  return (id, 1)
end
