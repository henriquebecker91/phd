module Kurpel2018

using JuMP

"""
    build(model, lwh′, XYZ, b, v)
    build(model, l′, w′, h′, X, Y, Z, b, v)

Add the variables and constraints of Kurpel 2018 3D-packing (minimize wasted
container volume, weakly heterogeneous containers, heterogeneous boxes,
free rotation allowed) to the model. The report presents formulations for
many variants using about the same main variables and constraints. The model
implemented here is not the same as any of the presented, but a slight
variation of the model in Section 2.1 (instead homogeneous containers and
minimizing the number of containers, here we have heterogeneous containers
and minimize the space wasted by the containers, or the total volume of the
containers used). This was done to allow comparison with the Chen 1995 model.
The parameter `b` (amount of boxes of some type) was removed, making the
formulation more BPP-like (instead of CSP-like).

For more information see: Kurpel, Deidson Vitorio, Cleder Marcos Schenekemberg,
Cassius Tadeu Scarpin, José Eduardo Pécora Junior, and Leandro C. Coelho.
“The Exact Solutions of Several Classes of Container Loading Problems,” n.d.

# Arguments
- `model`: JuMP model to which the variables and contraints will be added to.
- `lwh′`: three sequences of numbers (with the same length) denoting
  respectively the length/width/height of the boxes/cartons
- `XYZ`: three sequences of numbers (with the same length) denoting
  respectively the length/width/height of the containers
"""
function build(model, l′, w′, h′, X, Y, Z)
  build(model, (l′, w′, h′), (X, Y, Z))
end
function build(model, lwh′, XYZ)
  @assert length(lwh′) == 3
  @assert (isone ∘ length ∘ unique ∘ map)(length, lwh′)
  L, W, H = lwh′ # temporary names while the l/w/h is not yet created
  m = length(L)

  # the true l/h/w parameters take orientation in account
  l = Matrix(undef, m, 6)
  w = Matrix(undef, m, 6)
  h = Matrix(undef, m, 6)
  for i = 1:m
    l[i,1] = L; l[i,2] = L; l[i,3] = W; l[i,4] = W; l[i,5] = H; l[i,6] = H
    w[i,1] = W; w[i,2] = H; w[i,3] = L; w[i,4] = H; w[i,5] = L; w[i,6] = W
    h[i,1] = H; h[i,2] = W; h[i,3] = H; h[i,4] = L; h[i,5] = W; h[i,6] = L
  end

  @assert length(XYZ) == 3
  @assert (isone ∘ length ∘ unique ∘ map)(length, XYZ)
  X, Y, Z = XYZ
  x, y, z = max(X), max(Y), max(Z)
  C = lenght(X)

  v = map((a,b,c) -> a*b*c, X, Y, Z)

  @variables model begin
    x[1:n, 1:6, 1:C, 1:x, 1:y, 1:z], Bin
    e[1:C]
  end

  # Objective funtion Eq. (12)
  #@objective(model, Max, sum(
  #  x[i,g,j,x,y,z] * v[i] for i = 1:m, g = 1:6, j = 1:C, x = 1:(X-l[i,g]+1),
  #    y = 1:(Y-w[i,g]+1), z = 1:(Z-h[i,g]+1)
  #))
  # Objective funtion Eq. (7)
  @objective(model, Min, sum(e[j] * v[j] for j = 1:C))

  # Constraint (13) -- no overlap.
  for j = 1:C, s = 1:X[j], t = 1:Y[j], u = 1:Z[j]
    @constraint(model, sum(
      x[i, g, p, q, r] for i = 1:m, g = 1:6, p = max(1,s-l[i,g]+1):s, 
        q = max(1,q-w[i,g]+1):q, r = max(1,r-h[i,g]+1):r
    ) <= e[j])
  end

  # Constraint (14) -- every box needs to be packed.
  for i = 1:m
    @constraint(model, sum(
      u[i,g,j,p,q,r] for g = 1:6, x = 1:X[j], y = 1:Y[j], z = 1:Z[j]
    ) == 1)
  end

  model
end

end # module Kurpel2018

