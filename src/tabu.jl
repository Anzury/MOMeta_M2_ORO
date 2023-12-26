include("moves.jl")

# Tabu method : returns an improved solution w.r.t to sol, uses
# a parameter TL for the size of the memory (Tenure) and k 
# representing the maximum number of iterations without improvement
function tabu(ins::Instance, sol::Solution; TL=7, k=4, firstobj=true, move=swap, getMoves=getNonTabuSwapMoves)
  bestValue = firstobj ? getZ1(ins, sol.x, sol.y, sol.z) : getZ2(ins, sol.x)
  bestx = copy(sol.x)
  besty = copy(sol.y)
  bestz = copy(sol.z)
  num_iter = 0
  noimprovementcounter = 0
  TM = Array{Tuple, 1}(undef, TL)
  for i in 1:TL
    TM[i] = (-1, -1)
  end

  while noimprovementcounter < k
    found_improvement = false
    moves = getMoves(ins, sol, TM, TL, bestValue, firstobj=firstobj)

    # swap move
    new_sol = move(ins, sol, moves[2:end]...)
    if firstobj
      if new_sol.z1 < bestValue
        bestValue = new_sol.z1
        bestx = copy(new_sol.x)
        besty = copy(new_sol.y)
        bestz = copy(new_sol.z)
        found_improvement = true
      end
    else
      if new_sol.z2 < bestValue
        bestValue = new_sol.z2
        bestx = copy(new_sol.x)
        besty = copy(new_sol.y)
        bestz = copy(new_sol.z)
        found_improvement = true
      end
    end

    TM[num_iter%TL+1] = moves
    num_iter += 1
    if !found_improvement
      noimprovementcounter += 1
    end
  end

  z1 = firstobj ? bestValue : getZ1(ins, bestx, besty, bestz)
  z2 = firstobj ? getZ2(ins, bestx) : bestValue

  return Solution(bestx, besty, bestz, z1, z2)
end

# Improve a population of solutions using tabu search
function tabuImprovement(ins::Instance, pop::Vector{Solution}; TL=7, k=4)
  halfpop = Int(ceil(length(pop)/2))
  for i in 1:halfpop
    pop[i] = tabu(ins, pop[i], TL=TL, k=k, move=swap, getMoves=getNonTabuSwapMoves)
    pop[i] = tabu(ins, pop[i], TL=TL, k=k, move=shift, getMoves=getNonTabuShiftMoves)
    pop[i] = tabu(ins, pop[i], TL=TL, k=k, move=drop, getMoves=getNonTabuDropMoves)
    pop[i] = tabu(ins, pop[i], TL=TL, k=k, move=openLvl1, getMoves=getNonTabuOpenLvl1Moves)
  end

  for i in halfpop+1:length(pop)
    pop[i] = tabu(ins, pop[i], firstobj=false, TL=TL, k=k, move=swap, getMoves=getNonTabuSwapMoves)
    pop[i] = tabu(ins, pop[i], firstobj=false, TL=TL, k=k, move=shift, getMoves=getNonTabuShiftMoves)
    pop[i] = tabu(ins, pop[i], firstobj=false, TL=TL, k=k, move=drop, getMoves=getNonTabuDropMoves)
    pop[i] = tabu(ins, pop[i], firstobj=false, TL=TL, k=k, move=openLvl1, getMoves=getNonTabuOpenLvl1Moves)
  end

  return ins, pop
end
