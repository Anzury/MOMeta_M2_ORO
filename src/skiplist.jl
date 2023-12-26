using SkipLists
import Base

Tuple{Float64, Float64}(x::Float64, y::Float64) = (x, y)

const Point = Tuple{Float64, Float64}

dom(a::Point, b::Point) = a[1] <= b[1] && a[2] <= b[2] && a != b

# push! method for SkipList :
# if p is not in s and it is not dominated by any point in s, then we add p to s
# if p dominates a point in s, then we remove the dominated point from s and we
# add p to s
# if p is dominated by a point in s, then we do nothing
# if p is in s, then we do nothing
function Base.push!(s::SkipListSet{Point}, p::Point)
  if p in s
    return
  end
  todelete = []
  pisdominated = false
  for q in s
    if dom(p, q)
      push!(todelete, q)
    end
    if dom(q, p)
      pisdominated = true
    end
  end

  for q in todelete
    delete!(s, q)
  end

  if !pisdominated
    insert!(s, p)
  end
end

function initSkipList()
  return SkipListSet{Point}()
end
