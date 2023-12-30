using GeoJSON, DataFrames, Random

include("skiplist.jl")

mutable struct Instance
  fname::String
  nTerm::Int64
  nConc::Int64
  nLvl1::Int64
  nLvl2::Int64
  locationLvl1::Array{Float32,2}
  locationLvl2::Array{Float32,2}
  terminals::Array{Float32,2}
  B::Matrix{Int64}
  C::Matrix{Int64}
  S::Vector{Int64}
end

mutable struct Solution
  x::Array{Int64, 2}
  y::Array{Int64, 2}
  z::Array{Int64, 1}
  z1::Float64
  z2::Float64
end

mutable struct SolutionSet
  instance::Instance
  result_count::Int64
  YN::Vector{Point}
  solutions::Vector{Solution}
end

function getInstanceName(fname::String)
  return string(split(split(fname, '/')[end], '.')[1])
end

# get the cost matrix (distances) between locations in A and B
function getDistanceMatrix(A::Array{Float32,2}, B::Array{Float32,2})
  n = size(A, 1)
  m = size(B, 1)
  D = zeros(Int64, n, m)

  for i=1:n
    for j=1:m
      dist = (A[i,1]-B[j,1])^2 + (A[i,2]-B[j,2])^2
      D[i,j] = trunc(Int, dist^0.5)
    end
  end

  return D
end

# load instance : if path provided is a single txt file then use loadInstanceTxt,
# if it is a folder then get the first three geojson files in lexicographical
# order and use loadInstanceGeoJSON
function loadInstance(path::String)
  if isfile(path)
    return loadInstanceTxt(path)
  elseif isdir(path)
    files = readdir(path)
    files = sort(files)
    if endswith(files[1], ".geojson")
      path1 = joinpath(path, files[1])
      path2 = joinpath(path, files[2])
      path3 = joinpath(path, files[3])
      insname = getInstanceName(path)
      return loadInstanceGeoJSON(path1, path2, path3, insname)
    else
      error("Directory provided does not contain geojson files")
    end
  else
    error("Path provided is neither a file nor a directory")
  end
end

#=
# --------------------------------------
#           LOAD TXT INSTANCE
# --------------------------------------
# Input   : filename
# Output  : m, n, nLvl1, nLvl2, locationLvl1, locationLvl2, terminals, C, B
#
# nLvl1 the number of level 1 concentrators
# locationLvl1 the locations of the level 1 concentrators
# nLvl2 and locationLvl2 the equivalent of nLvl1 and locationLvl1 but for level
#   2 concentrators
# terminals the locations of the terminals
# C the cost matrix (distance matrix) between concentrators
# B the cost matrix (distance matrix) between terminals and concentrators
#           -----------------           
# File format :
#   - the first line contains four integers, which corresponds to m (# of
#   concentrators), n (# of terminals), p and r (both not used here)
#   - the next m lines represent the potential locations (x and y two floats
#   corresponding to the coordinates) of concentrators
#   - the next n lines respresent the locations (x and y two floats
#   corresponding to the coordinates) of the terminals and the weight
#   associated to that terminal (1 in all the examples considered here)
=#
function loadInstanceTxt(fname::String)
  f = open(fname)
  # read of the first line
  m::Int64, n::Int64, p::Int64, r::Int64 = parse.(Int, split(readline(f)))

  # we want to divide the set of concentrators in two levels so:
  # we assign 4/5 to lvl 1 concentrators and 1/5 to lvl 2 concentrators
  ratio::Float32 = 4/5
  nLvl1::Int64 = (Int)(m*ratio)
  nLvl2::Int64 = m - nLvl1

  # read the location of the lvl 1/2 concentrators and terminals
  locationLvl1 = stack([parse.(Float32, split(readline(f))) for i=1:nLvl1], dims=1)
  locationLvl2 = stack([parse.(Float32, split(readline(f))) for i=1:nLvl2], dims=1)
  terminals = stack([parse.(Float32, split(readline(f))) for i=1:n], dims=1)

  close(f)
  
  B = getDistanceMatrix(locationLvl1, locationLvl2)
  C = getDistanceMatrix(terminals, locationLvl1)
  mindist, maxdist = minimum(B), maximum(B)
  Random.seed!(42)
  S = rand(mindist:maxdist, nLvl2)
  Random.seed!()

  return Instance(fname, n, m, nLvl1, nLvl2, locationLvl1, locationLvl2, terminals, B, C, S)
end

#=
# --------------------------------------
#          LOAD GEOJSON INSTANCE
# --------------------------------------
# Input   : lvl1filename, lvl2filename, terminalsfilename
# Output  : m, n, nLvl1, nLvl2, locationLvl1, locationLvl2, terminals, C, B
#
# nLvl1 the number of level 1 concentrators
# locationLvl1 the locations of the level 1 concentrators
# nLvl2 and locationLvl2 the equivalent of nLvl1 and locationLvl1 but for level
#   2 concentrators
# terminals the locations of the terminals
# C the cost matrix (distance matrix) between concentrators
# B the cost matrix (distance matrix) between terminals and concentrators
#           -----------------           
# File format :
#   - GeoJSON files each containing the coordinates of the points
=#
function loadInstanceGeoJSON(lvl1fname::String, lvl2fname::String, termfname::String, insname::String="")
  function getCoordinates(geojson::DataFrame)
    # number of GeoJSON.Point in the GeoJSON instance
    nPoints = sum([g isa GeoJSON.Point for g in geojson[!, :geometry]])

    coords, i = zeros(Float32, nPoints, 2), 1
    for g in geojson[!, :geometry]
      if g isa GeoJSON.Point
        for j=1:2
          coords[i,j] = g.coordinates[j]
        end
        i += 1
      end
    end

    return coords
  end

  if insname == ""
    insname = "geojson"*randstring(6)
  end

  parseGeoJSON(filename) = DataFrame(GeoJSON.read(filename))

  locationLvl1 = getCoordinates(parseGeoJSON(lvl1fname))
  locationLvl2 = getCoordinates(parseGeoJSON(lvl2fname))
  terminals = getCoordinates(parseGeoJSON(termfname))

  nLvl1 = size(locationLvl1, 1)
  nLvl2 = size(locationLvl2, 1)
  n = size(terminals, 1)
  m = nLvl1 + nLvl2

  B = getDistanceMatrix(locationLvl1, locationLvl2)
  C = getDistanceMatrix(terminals, locationLvl1)
  mindist, maxdist = minimum(B), maximum(B)
  Random.seed!(42)
  S = rand(mindist:maxdist, nLvl2)
  Random.seed!()

  return Instance(insname, n, m, nLvl1, nLvl2, locationLvl1, locationLvl2, terminals, B, C, S)
end

# returns the value of a solution (x, y, z) for the 1st obj function
function getZ1(ins::Instance, x::Array{Int64, 2}, y::Array{Int64, 2}, z::Array{Int64, 1})
  return sum(ins.C[i,j]*x[i,j] for j in 1:ins.nLvl1 for i in 1:ins.nTerm) +
         sum(ins.B[j,k]*y[j,k] for k in 1:ins.nLvl2 for j in 1:ins.nLvl1) +
         sum(ins.S[k]*z[k] for k in 1:ins.nLvl2) + 0.0
end

# returns the value of a solution (x, y, z) for the 2nd obj function
function getZ2(ins::Instance, x::Array{Int64, 2})
  return maximum([x[i,j]*ins.C[i,j] for i in 1:ins.nTerm, j in 1:ins.nLvl1]) + 0.0
end

# utility function used to know if there is a non-null value in the vector
function isconnected(v::Array{Int64, 1})
  for i in v
    if i != 0
      return true
    end
  end
  return false
end
