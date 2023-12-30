function getIdealAndNadir(skiplist)
  idealz1, idealz2, nadirz1, nadirz2 = Inf, Inf, -1, -1
  for (z1, z2) in skiplist
    idealz1 = min(idealz1, z1)
    idealz2 = min(idealz2, z2)
    nadirz1 = max(nadirz1, z1)
    nadirz2 = max(nadirz2, z2)
  end

  return (idealz1, idealz2), (nadirz1, nadirz2)
end

function getLocalIdealAndNadir(p1, p2)
  return (min(p1[1], p2[1]), min(p1[2], p2[2])), (max(p1[1], p2[1]), max(p1[2], p2[2]))
end

function computeArea2points(p1, p2)
  return abs(p1[1] - p2[1]) * abs(p1[2] - p2[2])
end

function computeArea(YN, YPN)
  area = 0

  if length(YPN) == 1
    idealYN, _ = getIdealAndNadir(YN)
    return computeArea2points(idealYN, YPN[1])
  end

  for i=1:length(YPN)-1
    area += computeArea2points(YPN[i], YPN[i+1])
  end
  return area
end

function computeAreaExact(YN, YPN)
  under = computeArea(YN, YPN)
  idealYN, _ = getIdealAndNadir(YN)
  _, nadirYPN = getIdealAndNadir(YPN)
  square = computeArea2points(idealYN, nadirYPN)

  return square - under
end

function computeDistance2points(p1, p2)
  return sqrt((p1[1] - p2[1])^2 + (p1[2] - p2[2])^2)
end

function computeIdealsDistance(YN, YPN)
  idealYN, _ = getIdealAndNadir(YN)
  idealYPN, _ = getIdealAndNadir(YPN)

  return computeDistance2points(idealYN, idealYPN)
end
