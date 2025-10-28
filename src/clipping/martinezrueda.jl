# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
  MartinezRuedaClipping()

The Martinez-Rueda algorithm for clipping concave polygons.

## References

* Martínez, F., Rueda, A.J., Feito, F.R. 2009. [A new algorithm for computing Boolean operations on
  polygons](https://doi.org/10.1016/j.cag.2009.03.003)

### Notes
The algorithm works for both convex and concave polygons.
"""
struct MartinezRuedaClipping <: ClippingMethod end

function clip(poly::Polygon, other::Geometry, method::MartinezRuedaClipping)
  c = [clip(ring, boundary(other), method) for ring in rings(poly)]
  # r = [r for r in c if !isnothing(r)]
  # isempty(r) ? nothing : PolyArea(r)
end

function clip(ring::Ring, other::Ring, ::MartinezRuedaClipping)
  # get segments from both rings
  segs₁ = collect(segments(ring))
  segs₂ = collect(segments(other))

  # find all self-intersection points
  ps₁, idxs₁ = pairwiseintersect(segs₁)
  ps₂, idxs₂ = pairwiseintersect(segs₂)
  # update segments to include self-intersections
  if !isempty(ps₁)
    _updateselfintersections(segs₁, ps₁, idxs₁)
  end
  n = length(segs₁)

  if !isempty(ps₂)
    _updateselfintersections(segs₂, ps₂, idxs₂)
  end

  allsegs = Iterators.flatten((segs₁, segs₂))

  intersections, seginds = pairwiseintersect(allsegs)

  newsegs = Segment[]
  for (point, segidxs) in zip(intersections, seginds)
  end

  # collect all relevant points (vertices + intersections)
  P = eltype(intersections)
  points = P[]

  # if no points found, no intersection
  isempty(points) && return nothing

  unique!(points)
  length(points) < 3 && return nothing

  # create result ring
  Ring(points)
end

function _updateselfintersections(segs, ps, idxs)
  # split segments at intersection points
  counter = 0
  for (point, segidxs) in zip(ps, idxs)
    popat!(segs, segidxs[1] + counter - 1)
    idxs = [idx + counter for idx in segidxs]
    for idx in idxs
      s = segs[idx]
      # endpoints
      p₁, p₂ = vertices(s)
      # create new segments
      newseg₁ = Segment(p₁, point)
      newseg₂ = Segment(point, p₂)
      insert!(segs, idx, newseg₁)
      insert!(segs, idx + 1, newseg₂)
      counter += 1
    end
  end
end
