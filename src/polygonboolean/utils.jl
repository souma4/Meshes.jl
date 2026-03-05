# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# group intersections by (ringind, segind) using global segment indices
# returns Vector{Dict{Int,Vector{P}}} - one dict per ring
function _groupbyring(intersections, seginds, ringlengths)
  P = eltype(intersections)
  nrings = length(ringlengths)
  offsets = [0; cumsum(ringlengths)]

  # one dict per ring: segind -> points
  groups = [Dict{Int,Vector{P}}() for _ in 1:nrings]

  for (p, globalsegs) in zip(intersections, seginds), gind in globalsegs
    rind = searchsortedfirst(offsets, gind) - 1
    lind = gind - offsets[rind]
    push!(get!(groups[rind], lind, P[]), p)
  end

  groups
end

# splice/append intersection points into vertex vector for one ring
# insertions: Dict mapping segment index -> points
function _insertintoring!(v, n, insertions::Dict)
  isempty(insertions) && return v

  offset = 0
  for lind in sort!(collect(keys(insertions)))
    pts = insertions[lind]
    i = lind + offset

    pₛ = v[i]
    pₑ = lind < n ? v[i + 1] : v[1]

    # filter endpoints
    filter!(p -> !isapprox(p, pₛ) && !isapprox(p, pₑ), pts)
    isempty(pts) && continue

    # sort L-R, then order by segment direction
    sort!(pts)
    pₛ > pₑ && reverse!(pts)

    # splice or append
    if lind < n
      splice!(v, (i + 1):i, pts)
    else
      append!(v, pts)
    end

    offset += length(pts)
  end

  v
end

"""
    _insertintersections!(rings, intersections, seginds)

Insert intersection points into rings based on global segment indices from
`pairwiseintersect`.

"""
function _insertintersections!(rings, intersections, seginds)
  isempty(intersections) && return rings

  ringlengths = nvertices.(rings)
  groups = _groupbyring(intersections, seginds, ringlengths)

  for (i, (ring, group)) in enumerate(zip(rings, groups))
    isempty(group) && continue
    v = vertices(ring)
    _insertintoring!(v, ringlengths[i], group)
  end

  rings
end
