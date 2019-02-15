module CornerPoints

using Check3DPackings

"""
  compute_free_cps(xyzpqr, LWH)
  compute_free_cps(x, y, z, p, q, r, L, W, H)

# Arguments
- `xyzpqr`: a sextuple of collections, the first three are the position of the
  boxes in the x/y/z axis, the last three are the sizes of the boxes in th
  x/y/z axis.
- `LWH`: is a triple of the size of the container in the x/y/z axis.
# Returns
  - A vector of triples of triples. The first triple has the x
# Examples
  
"""
function compute_free_cps(xyzpqr, LWH)
  compute_free_cps(xyzpqr..., LWH...)
end
function compute_free_cps(x, y, z, p, q, r, L, W, H)
  xyzpqr = (x, y, z, p, q, r)
  @assert (isone ∘ length ∘ unique ∘ map)(length, xyzpqr)
  @assert all(map(a -> all(v -> v >= 1, a), xyzpqr))

  # Create the container walls.
  map(a->push!(a...), zip(xyzpqr, (0, 1, 1, 1, W, H)))
  map(a->push!(a...), zip(xyzpqr, (1, 0, 1, L, 1, H)))
  map(a->push!(a...), zip(xyzpqr, (1, 1, 0, L, W, 1)))
  n = length(x)
  
  # Change of notation for arrays with container walls.
  l = p; w = q; h = r

  # returned data structure
  xyzpqr′ijks = (Vector{typeof(L)}[], Vector{typeof(L)}[], Vector{typeof(L)}[],
    Vector{typeof(L)}[], Vector{typeof(L)}[], Vector{typeof(L)}[],
    Vector{Int64}[], Vector{Int64}[], Vector{Int64}[])

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
        p′ = max(x′, x[k], x[j])
        q′ = max(y′, y[k], y[i])
        r′ = max(z′, z[j], z[i])
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

  xyzpqr′ijk
end

end

