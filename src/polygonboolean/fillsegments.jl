# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    FillSegments{S}

Container for segments with fill annotations for geometric boolean algorithm.
"""
mutable struct FillSegments{S}
  segments::Vector{S}
  # fill status (subject/clip, filled above/below)
  subjectabove::BitVector
  subjectbelow::BitVector
  clipabove::BitVector
  clipbelow::BitVector
  # segment has been processed?
  computed::BitVector
end

function FillSegments(segments)
  n = length(segments)
  FillSegments(collect(segments), falses(n), falses(n), falses(n), falses(n), falses(n))
end

Base.collect(fs::FillSegments) = fs.segments
Base.length(fs::FillSegments) = length(fs.segments)
Base.getindex(fs::FillSegments, i::Int) = fs.segments[i]
Base.iterate(fs::FillSegments, state...) = iterate(fs.segments, state...)

setfills!(fs::FillSegments, above, below, i, issubject) =
  issubject ? _setsubjectfills!(fs, above, below, i) : _setclipfills!(fs, above, below, i)

_setsubjectfills!(fs::FillSegments, above, below, i) = (fs.subjectabove[i] = above; fs.subjectbelow[i] = below)

_setclipfills!(fs::FillSegments, above, below, i) = (fs.clipabove[i] = above; fs.clipbelow[i] = below)
getfills(fs::FillSegments, i, issubject) =
  issubject ? (fs.subjectabove[i], fs.subjectbelow[i]) : (fs.clipabove[i], fs.clipbelow[i])

_getsubjectfills(fs::FillSegments, i) = (fs.subjectabove[i], fs.subjectbelow[i])
_getclipfills(fs::FillSegments, i) = (fs.clipabove[i], fs.clipbelow[i])

_setcomputed!(fs::FillSegments, i::Int) = (fs.computed[i] = true)

"""
    segmentisless(fs, i, j)

Compare two segments by index for sweep line ordering.
Only works for segments that have been split at intersections.
"""
function segmentisless(fs::FillSegments, i::Int, j::Int)
  segi, segj = fs[i], fs[j]
  a1, a2 = sort(vertices(segi))
  b1, b2 = sort(vertices(segj))

  s1 = sideof(a1, Line(b1, b2))
  s2 = sideof(a2, Line(b1, b2))

  # RIGHT means a is below b
  if s1 == ON
    # if collinear, use endpoints, then index to break tie, otherwise use side of second point and b2 below a2
    s2 == ON ? (a1 < b1 || a2 < b2 || i > j) : s2 == RIGHT
  else
    # if both on same side, return that, otherwise use side of b
    s1 == s2 ? s1 == RIGHT : sideof(b1, Line(a1, a2)) != RIGHT
  end
end

# wrapper for AVL tree key sorting
struct SegmentIndex{S}
  fs::FillSegments{S}
  ind::Int
end

Base.isless(a::SegmentIndex, b::SegmentIndex) = segmentisless(a.fs, a.ind, b.ind)

# helper for storing event points
struct SegmentEvent
  ind::Int
  isstart::Bool
  issubject::Bool
end

# segment fill bitmasks
# describe whether segment is filled by subject or clip polygon above or below
const NONE = 0b0000
const SUBJTOP = 0b0001
const SUBJBOTTOM = 0b0010
const CLIPTOP = 0b0100
const CLIPBOTTOM = 0b1000
const BOTHTOP = SUBJTOP | CLIPTOP
const BOTHBOTTOM = SUBJBOTTOM | CLIPBOTTOM

"""
    _filled(operation, bits, subjmask, clipmask)

Determine if segment is filled based on operation and fill bits.

Fill rules:
- `union`: filled if either subject or clip is filled
- `intersect`: filled if both are filled
- `setdiff`: filled if subject filled and clip not
- `symdiff`/`xor`: filled if exactly one is filled
"""
function _filled(operation, bits, subjmask, clipmask)
  sj = (bits & subjmask) != 0
  cl = (bits & clipmask) != 0

  if operation == intersect
    sj && cl
  elseif operation == union
    sj || cl
  elseif operation == setdiff
    sj && !cl
  elseif operation == symdiff || operation == xor
    sj ⊻ cl
  else
    operation(sj, cl)
  end
end

# convenience wrappers for above/below fill checks
_filledabove(operation, bits) = _filled(operation, bits, SUBJTOP, CLIPTOP)
_filledbelow(operation, bits) = _filled(operation, bits, SUBJBOTTOM, CLIPBOTTOM)
