module CornerPoints

using Check3DPackings

export compute_free_cps

"""
  compute_free_cps(xyzpqr, LWH)
  compute_free_cps(x, y, z, p, q, r, L, W, H)

# Arguments
- `xyzpqr`: a sextuple of collections, the first three are the position of the
  boxes in the x/y/z axis, the last three are the sizes of the boxes in th
  x/y/z axis.
- `LWH`: is a triple of the size of the container in the x/y/z axis.
# Returns
- A nonuple of vectors. The first three vectors are the x/y/z positions of
  the corner points. The fourth, fifth, and sixth vectors are the minimum
  length/width/height a box needs to have to be placed in that corner point. 
  The last three vectors have the indexes of the boxes which provide the
  x/y/z-planes that form the corner point. The n+1, n+2, and n+3 values
  in those last three vectors denote, respectively, the back, left, and
  bottom faces of the container (that do not pertain to any box but are
  also used to form corner points).
# Examples
  
"""
function compute_free_cps(xyzpqr, LWH)
  compute_free_cps(xyzpqr..., LWH...)
end

function compute_free_cps(x, y, z, p, q, r, L, W, H)
  xyzpqr = (x, y, z, p, q, r)
  @assert (isone ∘ length ∘ unique ∘ map)(length, xyzpqr)
  @assert all(map(a -> all(v -> v >= 1, a), xyzpqr))
  for i = 1:length(x), j = (i+1):length(x)
    @assert !has_overlap(
      x[i], y[i], z[i], p[i], q[i], r[i],
      x[j], y[j], z[j], p[j], q[j], r[j]
    )
  end

  # Create the container walls.
  map(a->push!(a...), zip(xyzpqr, (0, 1, 1, 1, W, H)))
  map(a->push!(a...), zip(xyzpqr, (1, 0, 1, L, 1, H)))
  map(a->push!(a...), zip(xyzpqr, (1, 1, 0, L, W, 1)))
  n = length(x)
  
  # Change of notation for arrays with container walls.
  l = p; w = q; h = r

  # returned data structure
  xyzpqr′ijks = (Vector{typeof(L)}(), Vector{typeof(L)}(), Vector{typeof(L)}(),
    Vector{typeof(L)}(), Vector{typeof(L)}(), Vector{typeof(L)}(),
    Vector{Int64}(), Vector{Int64}(), Vector{Int64}())

  for k = 1:n
    z′ = z[k] + h[k]
    for j = 1:n
      y′ = y[j] + w[j]
      (z[j] + h[j] <= z′ || y[k] + w[k] <= y′) && continue
      for i = 1:n
        x′ = x[i] + l[i]
        ( x[k] + l[k] <= x′ ||
          x[j] + l[j] <= x′ ||
          z[i] + h[i] <= z′ ||
          y[i] + w[i] <= y′) && continue
        p′ = max(x′, x[k], x[j]) - x′ + 1
        q′ = max(y′, y[k], y[i]) - y′ + 1
        r′ = max(z′, z[j], z[i]) - z′ + 1
        for o = 1:(n-3) # the walls can never block a corner point
          has_overlap(
            x[o], y[o], z[o], l[o], w[o], h[o], x′, y′, z′, p′, q′, r′
          ) && @goto after_pushing_cp
        end
        xyzpqr′ijk = (x′, y′, z′, p′, q′, r′, i, j, k)
        map(a->push!(a...), zip(xyzpqr′ijks, xyzpqr′ijk))
        @label after_pushing_cp
      end
    end
  end

  # Remove the container walls from the list of boxes.
  map(a -> resize!(a, n - 3), xyzpqr)

  xyzpqr′ijks
end

end

