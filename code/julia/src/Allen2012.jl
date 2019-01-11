module Allen2012

using JuMP

# CODE ABANDONED: REASON: Seems like the family of constraints (17) that
# should prevent overlapping boxes does not work. There are more recent
# works with formulations that does not seem to share the same mistake,
# so we will use one of those.

"""
    build(model, pqr, LWH; allow_rot = false)
    build(model, p, q, r, L, W, H; allow_rot = false)

Add the variables and constraints of Allen 2012 3D-packing (maximize profit
of packed items, one single container, weakly heterogenous boxes, free rotation
can be allowed or not) to the model.

For more information see: Allen, Sam D., Edmund K. Burke, and Jakub Mareček.
“A Space-Indexed Formulation of Packing Boxes into a Larger Box.” Operations
Research Letters 40, no. 1 (January 1, 2012): 20–24.
https://doi.org/10.1016/j.orl.2011.10.008.

# Arguments
- `model`: JuMP model to which the variables and contraints will be added to.
- `pqr`: three sequences of numbers (with the same length) denoting
  respectively the length/width/height of the boxes/cartons
- `LWH`: three numbers denoting respectively the length/width/height of the
  container
- `A`: a sequence of numbers (with the same lenght as p/q/r) denoting the
  amount of each item available for packing
- `w`: a sequence of numbers (with the same lenght as p/q/r) denoting the
  worth/profit of each item
- `allow_rot`: if the model should use orientation variables
  TODO: add the name of the orientation variables here
"""
function build(model, p, q, r, L, W, H, A, w; allow_rot = false)
  build(model, (p, q, r), (L, W, H), A, w; allow_rot = allow_rot)
end
function build(model, pqr, LWH, A, w; allow_rot = false)
  # TODO: allow the caller to ask for rotation
  @assert !allow_rot

  @assert length(pqr) == 3
  @assert (isone ∘ length ∘ unique ∘ map)(length, pqr)
  L1, L2, L3 = pqr
  n = length(L1)
  @assert length(A) == N
  @assert length(w) == N
  #@assert length(LWH) == 3
  #@assert (isone ∘ length ∘ unique ∘ map)(length, LWH)
  DX, DY, DZ = LWH
  #m = length(LWH[1])
  #container_volumes = map(prod, zip(LWH...))
  # highest length/width/height of a container
  #ML, MW, MH = map((xs)->max(xs...), LWH)
  #ML, MW, MH = (10000, 10000, 10000)
  
  # To avoid creating unnecessary constraints, and also to simplify the code,
  # two lists are created. One contains the invalid positions and another the
  # valid ones. The objective function and constraints iterate such lists.
  invalid = Vector{NTuple{4,Int64}}()
  valid = Vector{NTuple{4,Int64}}()
  for t = 1:n, x = 1:DX, y = 1:DY, z = 1:DZ
    if (x > DX - L1[t] + 1) || (y > DY - L2[t] + 1) || (z > DZ - L3[t] + 1)
      push!(invalid, (t, x, y, z))
    else
      push!(valid, (t, x, y, z))
    end
  end
  
  @variables model begin
    u[1:n, 1:DX, 1:DY, 1:DZ], Bin
  end

  # Objective funtion Eq. (13)
  @objective(model, Max, sum(u[t,x,y,z] * w[t] for (t, x, y, z) = valid))

  # Constraint (14) -- the lower-left-bottom corner of two or more boxes do
  # not share the same point.
  for x = 1:DX, y = 1:DY, z = 1:DZ
    @constraint(model, sum(u[t, x, y, z] for t = 1:n) <= 1)
  end

  # Constraint (15) -- u[t,x,y,z] == 0 for all (x + L_{1t} > D_X) OR
  # (y + L_{2t} > D_Y) OR (z + L_{3t} > D_Z); i.e., the boxes stay inside
  # the container.
  for (t, x, y, z) = invalid
    setlowerbound(u[t, x, y, z], 0.0)
    setupperbound(u[t, x, y, z], 0.0)
  end

  # Constraint (16) -- we do not pack more boxes of some type than the ones
  # available.
  for t = 1:n
    @constraint(model,
      sum(u[t, x, y, z] for x = 1:DX, y = 1:DY, z = 1:DZ) <= A[t]
    )
  end

  # Constraint (17) -- boxes do not overlap.
  for (t, x, y, z) = valid
    @constraint(model, sum(u[t, x′, y′, z′] for x′ = x:(x + L1[t]),
      y′ = y:(y + L2[t]), z′ = z:(z + L2[t])) <= 1)
  end

  model
end

