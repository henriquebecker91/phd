module Chen1995

using JuMP
#using GLPKMathProgInterface # was used in solve

"""
    build(model, pqr, LWH; only_lxlzwyhz = true)
    build(model, p, q, r, L, W, H; only_lxlzwyhz = true)

Add the variables and constraints of Chen 1995 3D-packing (minimize wasted
container volume, heterogeneous containers, heterogeneous boxes, free rotation
allowed) to the model.

For more information see: Chen, C. S., S. M. Lee, and Q. S. Shen. "An Analytical Model for the Container Loading Problem." European Journal of Operational Research 80, no. 1 (January 5, 1995): 68–76. https://doi.org/10.1016/0377-2217(94)00002-T.

# Arguments
- `model`: JuMP model to which the variables and contraints will be added to.
- `pqr`: three sequences of numbers (with the same length) denoting
  respectively the length/width/height of the boxes/cartons
- `LWH`: three sequences of numbers (with the same length) denoting
  respectively the length/width/height of the containers
- `only_lxlzwyhz`: if the model should use only the orientation variables
  lx, lz, wy, and hz (fast but arcane), or all of them (legible but slow).
"""
function build(model, p, q, r, L, W, H; only_lxlzwyhz = true)
  build(model, (p, q, r), (L, W, H); only_lxlzwyhz = only_lxlzwyhz)
end
function build(model, pqr, LWH; only_lxlzwyhz = true)
  @assert length(pqr) == 3
  @assert (isone ∘ length ∘ unique ∘ map)(length, pqr)
  p, q, r = pqr
  N = length(pqr[1])
  @assert length(LWH) == 3
  @assert (isone ∘ length ∘ unique ∘ map)(length, LWH)
  L, W, H = LWH
  m = length(LWH[1])
  container_volumes = map(prod, zip(LWH...))
  # highest length/width/height of a container
  #ML, MW, MH = map((xs)->max(xs...), LWH)
  ML, MW, MH = (10000, 10000, 10000)
  
  @variables model begin
    s[1:N, 1:m], Bin
    n[1:m], Bin
    x[1:N] >= 0
    y[1:N] >= 0
    z[1:N] >= 0
    lx[1:N], Bin
    lz[1:N], Bin
    wy[1:N], Bin
    hz[1:N], Bin
    a[1:N, 1:N], Bin
    b[1:N, 1:N], Bin
    c[1:N, 1:N], Bin
    d[1:N, 1:N], Bin
    e[1:N, 1:N], Bin
    f[1:N, 1:N], Bin
  end
  if !only_lxlzwyhz
    @variables model begin
      ly[1:N], Bin
      wx[1:N], Bin
      wz[1:N], Bin
      hx[1:N], Bin
      hy[1:N], Bin
    end
  end

  # disable all unused variables (the model only access positions in
  # which i < k, so we can delete all i >= k)
  for var = [:a, :b, :c, :d, :e, :f], i = 1:N, k = 1:i
    setlowerbound(model[var][i, k], 0.0)
    setupperbound(model[var][i, k], 0.0)
  end

  @objective(model, Min, sum(container_volumes[j] * n[j] for j = 1:m))

  for k = 1:N, i = 1:(k-1) # forall i < k ∈ N
    # Constraints (1)-(6) -- no overlap among boxes inside the same container.
    if only_lxlzwyhz
      @constraints model begin
        x[i] + p[i]*lx[i] + q[i]*(lz[i] - wy[i] + hz[i]) + r[i]*(1 - lx[i] - lz[i] + wy[i] - hz[i]) <= x[k] + (1 - a[i, k])*ML
        # same as above but swapping [i] and [k], and replacing 'a' by 'b'
        x[k] + p[k]*lx[k] + q[k]*(lz[k] - wy[k] + hz[k]) + r[k]*(1 - lx[k] - lz[k] + wy[k] - hz[k]) <= x[i] + (1 - b[i, k])*ML

        y[i] + q[i]*wy[i] + p[i]*(1 - lx[i] - lz[i]) + r[i]*(lx[i] + lz[i] - wy[i]) <= y[k] + (1 - c[i, k])*MW
        # same as above but swapping [i] and [k], and replacing 'c' by 'd'
        y[k] + q[k]*wy[k] + p[k]*(1 - lx[k] - lz[k]) + r[k]*(lx[k] + lz[k] - wy[k]) <= y[i] + (1 - d[i, k])*MW

        z[i] + r[i]*hz[i] + q[i]*(1 - lz[i] - hz[i]) + p[i]*lz[i] <= z[k] + (1 - e[i, k])*MH
        # same as above but swapping [i] and [k], and replacing 'e' by 'f'
        z[k] + r[k]*hz[k] + q[k]*(1 - lz[k] - hz[k]) + p[k]*lz[k] <= z[i] + (1 - f[i, k])*MH
      end
    else
      @constraints model begin
        x[i] + p[i]*lx[i] + q[i]*wx[i] + r[i]*hx[i] <= x[k] + (1 - a[i, k])*ML
        x[k] + p[k]*lx[k] + q[k]*wx[k] + r[k]*hx[k] <= x[i] + (1 - b[i, k])*ML
        y[i] + q[i]*wy[i] + p[i]*ly[i] + r[i]*hy[i] <= y[k] + (1 - c[i, k])*MW
        y[k] + q[k]*wy[k] + p[k]*ly[k] + r[k]*hy[k] <= y[i] + (1 - d[i, k])*MW
        z[i] + r[i]*hz[i] + q[i]*wz[i] + p[i]*lz[i] <= z[k] + (1 - e[i, k])*MH
        z[k] + r[k]*hz[k] + q[k]*wz[k] + p[k]*lz[k] <= z[i] + (1 - f[i, k])*MH
      end
    end

    for j = 1:m # forall i < k ∈ N, j ∈ m
      # Constraint (7) -- boxes in the same container must have a relative
      # position to each other.
      @constraint(model, 
        a[i, k] + b[i, k] + c[i, k] + d[i, k] + e[i, k] + f[i, k] >= s[i, j] + s[k, j] - 1
      )
    end
  end # end forall i < k ∈ N

  if !only_lxlzwyhz
    for i = 1:N
      # Necessary if the model includes all orientation variables.
      @constraints model begin
        lx[i] + ly[i] + lz[i] == 1
        wx[i] + wy[i] + wz[i] == 1
        hx[i] + hy[i] + hz[i] == 1
        lx[i] + wx[i] + hx[i] == 1
        ly[i] + wy[i] + hy[i] == 1
        lz[i] + wz[i] + hz[i] == 1
      end
    end
  end

  for i = 1:N
    # Constraint (8) -- every box must be inside one container.
    @constraint(model, sum(s[i, j] for j = 1:m) == 1)
  end

  # Constraint (9) -- if there is a box inside a container, the container
  # must be marked as used.
  sorted_bvs = sort!(map((xs)->*(xs...), zip(pqr...)))
  for j = 1:m
    # MJ is the basic volumetric upper bound in the amount of boxes that can be
    # fitted inside the respective container.
    #MJ = 0
    MJ = 10000
    #remaining_cv = container_volumes[j]
    #while MJ < N && remaining_cv >= sorted_bvs[MJ+1]
    #  remaining_cv -= sorted_bvs[MJ+1]
    #  MJ += 1
    #end
    @constraint(model, sum(s[i, j] for i = 1:N) <= MJ*n[j])
  end

  # Constraints (10)-(12) -- the boxes must respect the container boundaries
  # (the back-left-bottom boundary is already enforced by the '>= 0' trivial
  # constraint).
  for i = 1:N, j = 1:m
    if only_lxlzwyhz
      @constraints model begin
        x[i] + p[i]*lx[i] + q[i]*(lz[i] - wy[i] + hz[i]) + r[i]*(1 - lx[i] - lz[i] + wy[i] - hz[i]) <= L[j] + (1 - s[i, j])*ML
        y[i] + q[i]*wy[i] + p[i]*(1 - lx[i] - lz[i]) + r[i]*(lx[i] + lz[i] - wy[i]) <= W[j] + (1 - s[i, j])*MW
        z[i] + r[i]*hz[i] + q[i]*(1 - lz[i] - hz[i]) + p[i]*lz[i] <= H[j] + (1 - s[i, j])*MH
      end
    else
      @constraints model begin
        x[i] + p[i]*lx[i] + q[i]*wx[i] + r[i]*hx[i] <= L[j] + (1 - s[i, j])*ML
        y[i] + q[i]*wy[i] + p[i]*ly[i] + r[i]*hy[i] <= W[j] + (1 - s[i, j])*MW
        z[i] + r[i]*hz[i] + q[i]*wz[i] + p[i]*lz[i] <= H[j] + (1 - s[i, j])*MH
      end
    end
  end

  model
end

"""
    extract_attributions(model) :: Dict{Int64, Vector{Int64}}

Returns a dictionary with the index of the containers as keys and the index
of the boxes inside that container as values.
"""
function extract_attributions(model) :: Dict{Int64, Vector{Int64}}
  s = getvalue(model[:s])
  N = size(s, 1)
  m = size(s, 2)
  d = Dict{Int64, Vector{Int64}}()
  # builds the dictionary from container index to box indexes
  for i = 1:N, j = 1:m
    if isone(round(s[i, j]))
      if haskey(d, j)
        push!(d[j], i)
      else
        d[j] = [i]
      end
    end
  end
  
  d
end

"""
    extract_positions(model) :: NTuple{3, Vector{Float64}}

Returns three arrays (x, y, z) of coordinates where the boxes were placed.
This model has the coordinates of he boxes as continuous variables. It is
tempting to truncate the values, as this would keep the validity of the
solution and make it more compatible with other models, but as the model
can be modified to worry about the gravity center, it is best kept as it is.
"""
function extract_positions(model) :: NTuple{3, Vector{Float64}}
  map(s ->getvalue(model[s]), (:x, :y, :z))
end

"""
    extract_oriented_sizes(model, pqr; only_lxlzwyhz = true) :: NTuple{3, Vector{Int64}}

Returns the oriented sizes of the boxes in the solved model.
"""
function extract_oriented_sizes(
  model, pqr; only_lxlzwyhz = true
) :: NTuple{3, Vector{Int64}}
  # computes the oriented pqr values (sizes of the boxes when oriented)
  N = length(model[:x])
  zeroes = repeat([zero(typeof(pqr[1][1]))], N)
  p′, q′, r′ = zeroes, copy(zeroes), copy(zeroes)
  p, q, r = pqr
  lx = getvalue(model[:lx])
  lz = getvalue(model[:lz])
  wy = getvalue(model[:wy])
  hz = getvalue(model[:hz])
  if !only_lxlzwyhz
    ly = getvalue(model[:ly])
    wx = getvalue(model[:wx])
    wz = getvalue(model[:wz])
    hx = getvalue(model[:hx])
    hy = getvalue(model[:hy])
  end
  if only_lxlzwyhz
    for i = 1:N
      p′[i] = p[i]*lx[i] + q[i]*(lz[i] - wy[i] + hz[i]) +
        r[i]*(1 - lx[i] - lz[i] + wy[i] - hz[i])
      q′[i] = p[i]*(1 - lx[i] - lz[i]) + q[i]*wy[i] +
        r[i]*(lx[i] + lz[i] - wy[i])
      r′[i] = p[i]*(lz[i]) + q[i]*(1 - lz[i] - hz[i]) + r[i]*hz[i]
    end
  else
    for i = 1:N
      p′[i] = p[i]*lx[i] + q[i]*wx[i] + r[i]*hx[i]
      q′[i] = p[i]*ly[i] + q[i]*wy[i] + r[i]*hy[i]
      r′[i] = p[i]*lz[i] + q[i]*wz[i] + r[i]*hz[i]
    end
  end

  (p′, q′, r′)
end

"""
    extract_solution(model, pqr; only_lxlzwyhz = true)

Returns the triple of the results of extract\\_{attributions, positions,
oriented\\_sizes} over the solved model.

The value of only_lxlzwyhz needs to be the same passed to Chen1995.build.
"""
function extract_solution(model, pqr; only_lxlzwyhz = true)
  d = extract_attributions(model)
  c = extract_positions(model)
  pqr′ = extract_oriented_sizes(model, pqr; only_lxlzwyhz = only_lxlzwyhz)

  (d, c, pqr′)
end

#= PROBLEM: should we treat the possibility of the model being unfeasible?
"""
    solve(pqr, LWH, solver = GLPKSolverMIP())

Solves the problem defined by pqr and LWH using the Chen 1995 model in the
specified solver, and returns a triple of extract_{attributions, positions,
oriented_sizes} (which is the extracted solution).
"""
function solve(pqr, LWH, solver = GLPKSolverMIP())
  model = Model(solver = solver)
  build(model, pqr, LWH)
  status = solve(model, suppress_warnings = true)
  if 
  extract_solution(model, pqr)
end=#

#=function chen1995_1to7(model, pqr, LWH)
  
end=#

end # module Chen1995

