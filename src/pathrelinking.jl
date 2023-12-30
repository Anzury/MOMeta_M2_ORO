include("parser.jl")

function pathrelinking(ins::Instance, pop::Vector{Solution})
  halfpop = Int(ceil(length(pop)/2))
  S1 = pop[1:halfpop]
  S2 = pop[halfpop:end]

  # Sort pop by objective value
  sort!(S1, by = x -> x.z1)
  sort!(S2, by = x -> x.z2)

  for i=1:10
    push!(S1, pop[rand(1:halfpop)])
    push!(S2, pop[rand(halfpop:end)])
  end

  solutions = []
  YPN = initSkipList()

  for (s1, s2) in zip(S1, S2)
    # s1 is the solution dominating the other one
    if s2.z1 < s1.z1 && s2.z2 < s1.z2
      s1, s2 = s2, s1
    end
    push!(YPN, (s1.z1, s1.z2))
    push!(YPN, (s2.z1, s2.z2))

    s1lvl1, s2lvl1 = [], []
    for j=1:ins.nLvl1
      for i=1:ins.nTerm
        if s1.x[i, j] == 1
          push!(s1lvl1, j)
        end
        if s2.x[i, j] == 1
          push!(s2lvl1, j)
        end
      end
    end

    In = setdiff(s2lvl1, s1lvl1)
    Out = setdiff(s1lvl1, s2lvl1)
    lIn, lOut = length(In), length(Out)

    nx, ny, nz = copy(s1.x), copy(s1.y), copy(s1.z)
    px, py, pz = copy(nx), copy(ny), copy(nz)
    pz1, pz2 = s1.z1, s1.z2

    while lIn > 0 && lOut > 0
      jo = pop!(Out)
      ji = pop!(In)
      lIn -= 1
      lOut -= 1
      newly_connected = []
      ji_isinserted = false

      # reconnect all terminals to their closest lvl1 concentrator but prevent
      # any terminal to connect to jo
      for i=1:ins.nTerm
        if nx[i, jo] == 1
          # find the closest lvl1 concentrator
          minDist, target_j = Inf, -1
          for j=1:ins.nLvl1
            if j != jo && nx[i, j] == 0
              if ins.C[i, j] < minDist
                minDist = ins.C[i, j]
                target_j = j
              end
            end
          end
          nx[i, jo] = 0
          nx[i, target_j] = 1
          if !ji_isinserted && target_j == ji
            ji_isinserted = true
          end
          push!(newly_connected, target_j)

          # disconnect lvl2 concentrator that were connected to jo
          for k=1:ins.nLvl2
            if ny[jo, k] == 1
              ny[jo, k] = 0
              # if the lvl2 concentrator is not connected to any lvl1 concentrator
              # then close it
              if !isconnected(ny[:, k])
                nz[k] = 0
              end
            end
          end
        end
      end

      if !ji_isinserted
        push!(newly_connected, ji)
      end

      # check if lvl1 concentrators are connected to their closest lvl2 concentrator
      # if not, connect them
      for j in newly_connected
        alreadyconnected = false
        minDist, target_k = Inf, -1
        for k=1:ins.nLvl2
          if ny[j, k] == 0
            if ins.B[j, k] < minDist
              minDist = ins.B[j, k]
              target_k = k
            end
          else
            alreadyconnected = true
            break
          end
        end
        if !alreadyconnected
          ny[j, target_k] = 1
          nz[target_k] = 1
        end
      end

      # add the new solution to the archive (YPN)
      nz1 = getZ1(ins, nx, ny, nz)
      nz2 = getZ2(ins, nx)
      
      if nz1 < pz1 && nz2 < pz2
        pz1, pz2 = nz1, nz2
        px, py, pz = copy(nx), copy(ny), copy(nz)

        push!(solutions, Solution(nx, ny, nz, nz1, nz2))
        push!(YPN, (nz1, nz2))
      else
        nx, ny, nz = copy(px), copy(py), copy(pz)
        nz1, nz2 = pz1, pz2
      end

    end
  end

  return YPN, solutions
end
