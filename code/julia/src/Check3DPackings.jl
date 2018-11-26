module Check3DPackings

export has_violations, all_violations, has_overlap, is_inside,
  is_orientation_of

"""
    has_violations(c, xyzpqr, LWH) :: Bool
    has_violations(c, x, y, z, p, q, r, L, W, H) :: Bool

Alias for !isempty(all_violations(c, xyzpqr, LWH)).

```jldoctest
julia> xyzpqr = ([0, 3], [0, 4], [0, 5], [4, 7], [5, 7], [6, 7]);

julia> LWH = ([7, 6], [7, 6], [7, 6]);

julia> has_violations((Dict(1 => [1, 2])), xyzpqr, LWH)
true

julia> has_violations((Dict(2 => [2])), xyzpqr, LWH)
true

julia> has_violations((Dict(1 => [1])), xyzpqr, LWH)
false
```
"""
function has_violations(c, xyzpqr, LWH) :: Bool
  !isempty(all_violations(c, xyzpqr, LWH))
end
function has_violations(c, x, y, z, p, q, r, L, W, H) :: Bool
  has_violations(c, (x, y, z, p, q, r), (L, W, H))
end

"""
    all_violations(c, xyzpqr, LWH) :: Vector{Tuple{Symbol,Any}}
    all_violations(c, x, y, z, p, q, r, L, W, H) :: Vector{Tuple{Symbol,Any}}

Checks the packings for violation of the contrainsts that are common to all
variants of 3D-packing: boxes must be completely inside its respective
containers; boxes inside the same container must not overlap.

The function does not check if all boxes are assigned to a container (as
knapsack variants often cannot pack all boxes). The function assumes the boxes
are already rotated to the orientation used in the solution.

# Arguments
- `c`: a collection in which the keys are container indexes and the values are
  collections of indexes of all the boxes inside the respective container.
- `xyzpqr`: a sextuple of collections, the first three are the position of the
  boxes in the x/y/z axis, the last three are the sizes of the boxes in th
  x/y/z axis.
- `LWH`: a triple of collections consisting in the size of the containers in
  the x/y/z axis.

# Returns
- `[(:outside, (j, i))]`: if box i is not completely inside container j.
- `[(:overlap, (j, i, k))]`: if box i and k overlap (inside container j).
- `[]`: if none of common constraints are violated.

# Examples
```jldoctest
julia> xyzpqr = ([0, 3], [0, 4], [0, 5], [4, 7], [5, 7], [6, 7]);

julia> LWH = ([7, 6], [7, 6], [7, 6]);

julia> all_violations((Dict(1 => [1, 2])), xyzpqr, LWH)
2-element Array{Tuple{Symbol,Any},1}:
 (:overlap, (1, 1, 2))
 (:outside, (1, 2))

julia> all_violations((Dict(2 => [2])), xyzpqr, LWH)
1-element Array{Tuple{Symbol,Any},1}:
 (:outside, (2, 2))

julia> all_violations((Dict(1 => [1])), xyzpqr, LWH)
0-element Array{Tuple{Symbol,Any},1}
```
"""
function all_violations(
  c, x, y, z, p, q, r, L, W, H
) :: Vector{Tuple{Symbol,Any}}
  all_violations(c, (x, y, z, p, q, r), (L, W, H))
end
function all_violations(c, xyzpqr, LWH) :: Vector{Tuple{Symbol,Any}}
  @assert length(xyzpqr) == 6
  @assert length(LWH) == 3
  @assert (isone ∘ length ∘ unique ∘ map)(length, xyzpqr)
  @assert (isone ∘ length ∘ unique ∘ map)(length, LWH)
  violations = []

  index_all(a, i) = map(b -> b[i], a)
  for (j, box_idxs) in c
    container = index_all(LWH, j)
    for i in 1:length(box_idxs)
      box1 = index_all(xyzpqr, box_idxs[i])
      if !is_inside(box1..., container...)
        push!(violations, (:outside, (j, box_idxs[i])))
      end
      for k in (i+1):length(box_idxs)
        box2 = index_all(xyzpqr, box_idxs[k])
        if has_overlap(box1..., box2...)
          push!(violations, (:overlap, (j, box_idxs[i], box_idxs[k])))
        end
      end
    end
  end

  violations
end

"""
    has_overlap(xyzpqr, xyzpqr′) :: Bool
    has_overlap(x, y, z, p, q, r, x′, y′, z′, p′, q′, r′) :: Bool

Returns true if the two positioned boxes overlap; returns false otherwise.

If the two boxes share a point, line, or plane, they are not considered
overlapping. They only overlap if a new box could be defined from the
intersection of both boxes (i.e., the two boxes share some volume).

# Examples
```jldoctest
julia> has_overlap(0, 0, 0, 4, 5, 6, 4, 5, 6, 10, 10, 10)
false

julia> has_overlap((0, 0, 0, 4, 5, 6), (3, 4, 5, 10, 10, 10))
true
```
"""
function has_overlap(x, y, z, p, q, r, x′, y′, z′, p′, q′, r′) :: Bool
  intersect((c, s, c′, s′)) = c < (c′ + s′) && (c + s) > c′
  all(intersect.([(x, p, x′, p′), (y, q, y′, p′), (z, r, z′, p′)]))
end

function has_overlap(xyzpqr, xyzpqr′) :: Bool
  x, y, z, p, q, r = xyzpqr
  x′, y′, z′, p′, q′, r′ = xyzpqr′
  has_overlap(x, y, z, p, q, r, x′, y′, z′, p′, q′, r′)
end


"""
    is_inside(xyzpqr, xyzpqr′) :: Bool
    is_inside(x, y, z, p, q, r, x′, y′, z′, p′, q′, r′) :: Bool
    is_inside(xyzpqr, LWH) :: Bool
    is_inside(x, y, z, p, q, r, L, W, H) :: Bool

Returns true if the first positioned box is completely inside the second
positioned box; returns false otherwise.

If the positions of the second box are ommited (i.e., xyzpqr′ is replaced by
favor of LWH) the origin (0, 0, 0) is assumed for the positions, and LWH is
used for the sizes.

# Examples
```jldoctest
julia> is_inside(1, 1, 1, 9, 9, 9, 0, 0, 0, 10, 10, 10)
true

julia> is_inside((0, 0, 0, 10, 10, 10), (0, 0, 0, 10, 10, 10))
true

julia> is_inside(1, 1, 1, 4, 5, 6, 4, 5, 5)
false

julia> is_inside((1, 1, 1, 4, 5, 6), (4, 5, 5))
false
```
"""
function is_inside(x, y, z, p, q, r, x′, y′, z′, p′, q′, r′) :: Bool
  (x >= x′ && y >= y′ && z >= z′ && (x + p) <= (x′ + p′) &&
    (y + q) <= (y′ + q′) && (z + r) <= (z′ + r′))
end

function is_inside(
  xyzpqr :: Tuple{Any, Any, Any, Any, Any, Any},
  xyzpqr′ :: Tuple{Any, Any, Any, Any, Any, Any}
) :: Bool
  x, y, z, p, q, r = xyzpqr
  x′, y′, z′, p′, q′, r′ = xyzpqr′
  is_inside(x, y, z, p, q, r, x′, y′, z′, p′, q′, r′)
end

function is_inside(x, y, z, p, q, r, L, W, H) :: Bool
  is_inside(x, y, z, p, q, r, 0, 0, 0, L, W, H)
end

function is_inside(
  xyzpqr :: Tuple{Any, Any, Any, Any, Any, Any},
  LWH :: Tuple{Any, Any, Any}
) :: Bool
  x, y, z, p, q, r = xyzpqr
  L, W, H = LWH
  is_inside(x, y, z, p, q, r, 0, 0, 0, L, W, H)
end

"""
    is_orientation_of(pqr, pqr′) :: Bool
    is_orientation_of(p, q, r, p′, q′, r′) :: Bool

Returns true if pqr is a permutation of the values of pqr′; returns false
otherwise.

The function expects two boxes, not an arbitrary number of boxes. Note that
is_orientation_of.(p, q, r, p′, q′, r′) gives the expected answer if every
parameter is an array, while is_orientation_of(pqr, pqr′) gives the expected
answer if both parameters are arrays of triples.

# Examples
```jldoctest
julia> is_orientation_of((1, 2, 3), (1, 2, 3))
true

julia> is_orientation_of((1, 1, 2), (1, 2, 2))
false

julia> is_orientation_of(1, 2, 3, 1, 3, 2)
true

julia> is_orientation_of(1, 2, 1, 2, 1, 1)
true

julia> is_orientation_of(1, 2, 3, 3, 2, 2)
false
```
"""
function is_orientation_of(pqr, pqr′) :: Bool
  is_orientation_of(pqr..., pqr′...)
end
function is_orientation_of(p, q, r, p′, q′, r′) :: Bool
  (p == p′ && q == q′ && r == r′) ||
  (p == p′ && q == r′ && r == q′) ||
  (p == q′ && q == p′ && r == r′) ||
  (p == r′ && q == q′ && r == p′) ||
  (p == q′ && q == r′ && r == p′) ||
  (p == r′ && q == p′ && r == q′)
end

end # module Check3DPackings

