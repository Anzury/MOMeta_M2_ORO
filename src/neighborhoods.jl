include("parser.jl")

# Swap move : swap the terminals assigned to 2 lvl1 concentrators
function swap(ins::Instance, sol::Solution, j1::Int64, j2::Int64)
  nx::Array{Int64, 2} = copy(sol.x)
  ny::Array{Int64, 2} = copy(sol.y)
  nz::Array{Int64, 1} = copy(sol.z)
  tempx::Array{Int64, 2} = zeros(Int64, ins.nTerm, ins.nLvl1)

  # copy temporarily the terminals assigned to j1
  for i in 1:ins.nTerm
    if nx[i, j1] == 1
      tempx[i, j1] = 1
    end
  end

  # assign to j1 the terminals assigned to j2
  for i in 1:ins.nTerm
    nx[i, j1] = nx[i, j2]
  end

  # assign to j2 the terminals that were assigned to j1 (in tempx)
  for i in 1:ins.nTerm
    nx[i, j2] = tempx[i, j1]
  end

  return Solution(nx, ny, nz, getZ1(ins, nx, ny, nz), getZ2(ins, nx))
end

# Shift move : assign one terminal i to another (opened or not) lvl1 concentrator
function shift(ins::Instance, sol::Solution, i::Int64, start_j::Int64, target_j::Int64)
  nx::Array{Int64, 2} = copy(sol.x)
  ny::Array{Int64, 2} = copy(sol.y)
  nz::Array{Int64, 1} = copy(sol.z)

  # if the terminal i assigned to concentrator start_j was the last one assigned
  # to it then we close the lvl2 concentrators connected to start_j
  if sum(nx[:, start_j]) == 1
    for k in 1:ins.nLvl2
      ny[start_j, k] = 0
      if !isconnected(ny[:, k])
        nz[k] = 0
      end
    end
  end

  # is the target concentrator already connected to lvl2 concentrator ?
  target_opened = isconnected(ny[target_j, :])

  # make shift
  nx[i, start_j] = 0
  nx[i, target_j] = 1

  # if the target concentrator wasn't connected to lvl2 concentrator then
  # we connect it to the closest lvl2 concentrator (check with ins.B[j,k])
  if !target_opened
    # find the closest lvl2 concentrator
    d, target_k = typemax(Int64), -1
    for k in 1:ins.nLvl2
      if ins.B[target_j, k] < d
        d = ins.B[target_j, k]
        target_k = k
      end
    end
    # connect the target concentrator to the closest lvl2 concentrator
    ny[target_j, target_k] = 1
    nz[target_k] = 1
  end

  return Solution(nx, ny, nz, getZ1(ins, nx, ny, nz), getZ2(ins, nx))
end

# Drop move : drop the terminals assigned to a lvl1 concentrator and
# connect each terminal to the closest lvl1 concentrator
function drop(ins::Instance, sol::Solution, j1::Int64)
  nx::Array{Int64, 2} = copy(sol.x)
  ny::Array{Int64, 2} = copy(sol.y)
  nz::Array{Int64, 1} = copy(sol.z)

  dropped = []

  # drop the terminals assigned to j1
  for i in 1:ins.nTerm
    if nx[i, j1] == 1
      nx[i, j1] = 0
      push!(dropped, i)
    end
  end

  # check if there was a lvl2 connected to j1
  if isconnected(ny[j1, :])
    # if so then we close the lvl2 concentrators connected to j1
    for k in 1:ins.nLvl2
      ny[j1, k] = 0
      if !isconnected(ny[:, k])
        nz[k] = 0
      end
    end
  end

  # find the closest non assigned lvl1 concentrator for each dropped terminal
  for i in dropped
    d, target_j = typemax(Int64), -1
    for j in 1:ins.nLvl1
      if nx[i, j] == 0 && ins.C[i, j] < d
        d = ins.C[i, j]
        target_j = j
      end
    end
    nx[i, target_j] = 1

    # if lvl1 connector is not connected to lvl2 concentrator then we connect
    # it to the closest lvl2 concentrator
    if !isconnected(ny[target_j, :])
      d, target_k = typemax(Int64), -1
      for k in 1:ins.nLvl2
        if ins.B[target_j, k] < d
          d = ins.B[target_j, k]
          target_k = k
        end
      end
      ny[target_j, target_k] = 1
      nz[target_k] = 1
    end
  end

  return Solution(nx, ny, nz, getZ1(ins, nx, ny, nz), getZ2(ins, nx))
end

# openLvl1 move : open a lvl1 concentrator, connect it to the closest lvl2 and
# if the terminals assigned to other lvl1 concentrators are closer to this
# lvl1 concentrator then disconnect them from their farthest lvl1 concentrator
# and connect them to this lvl1 concentrator
function openLvl1(ins::Instance, sol::Solution, j1::Int64)
  nx::Array{Int64, 2} = copy(sol.x)
  ny::Array{Int64, 2} = copy(sol.y)
  nz::Array{Int64, 1} = copy(sol.z)

  candidate_terms = []
  terms = [i for i in 1:ins.nTerm if isconnected(nx[i, :])]
  # find all terminals that are closer to j1 than to their farthest assigned
  # lvl1 concentrator
  for i in terms
    first_j, second_j = -1, -1
    # get index of the two lvl1 concentrators connected to terminal i
    for j in 1:ins.nLvl1
      if nx[i, j] == 1 && j != j1
        if first_j == -1
          first_j = j
        else
          second_j = j
        end

        if second_j != -1
          break
        end
      end
    end
    if first_j == -1 || second_j == -1
      continue
    end
    farthest = ins.C[i, first_j] > ins.C[i, second_j] ? first_j : second_j

    # check if terminal i is closer to j than to its farthest assigned lvl1
    # concentrator
    if ins.C[i, j1] < ins.C[i, farthest]
      push!(candidate_terms, (i, farthest))
    end
  end

  for (i, j) in candidate_terms
    # disconnect terminal i from its farthest assigned lvl1 concentrator
    nx[i, j] = 0
    #check if terminal i was the last terminal assigned to j, if so then
    #close the lvl2 concentrators connected to j
    if !isconnected(nx[:, j])
      for k in 1:ins.nLvl2
        ny[j, k] = 0
        if !isconnected(ny[:, k])
          nz[k] = 0
        end
      end
    end
    nx[i, j1] = 1 # connect terminal i to j1
  end

  # connect j1 to the closest lvl2 concentrator
  if length(candidate_terms) > 0
    d, target_k = typemax(Int64), -1
    for k in 1:ins.nLvl2
      if ins.B[j1, k] < d
        d = ins.B[j1, k]
        target_k = k
      end
    end
    ny[j1, target_k] = 1
    nz[target_k] = 1
  end

  return Solution(nx, ny, nz, getZ1(ins, nx, ny, nz), getZ2(ins, nx))
end
