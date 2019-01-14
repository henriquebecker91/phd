module Kurpel2018

using JuMP

"""
    build(model, LWH, XYZ)
    build(model, L, W, H, X, Y, Z)

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
function build(model, L, W, H, X, Y, Z)
  build(model, (L, W, H), (X, Y, Z))
end
function build(model, LWH, XYZ)
  @assert length(LWH) == 3
  @assert (isone ∘ length ∘ unique ∘ map)(length, LWH)
  L, W, H = LWH # temporary names while the l/w/h is not yet created
  m = length(L)

  # the true l/h/w parameters take orientation in account
  l = Matrix(undef, m, 6)
  w = Matrix(undef, m, 6)
  h = Matrix(undef, m, 6)
  for i = 1:m
    a, b, c = L[i], W[i], H[i]
    l[i,1] = a; l[i,2] = a; l[i,3] = b; l[i,4] = b; l[i,5] = c; l[i,6] = c
    w[i,1] = b; w[i,2] = c; w[i,3] = a; w[i,4] = c; w[i,5] = a; w[i,6] = b
    h[i,1] = c; h[i,2] = b; h[i,3] = c; h[i,4] = a; h[i,5] = b; h[i,6] = a
  end

  @assert length(XYZ) == 3
  @assert (isone ∘ length ∘ unique ∘ map)(length, XYZ)
  X, Y, Z = XYZ

  C = length(X)

  v = map((a,b,c) -> a*b*c, X, Y, Z)

  @variables model begin
    # TODO: for each container dimension, x uses the largest value of that
    # dimension among all containers. This way, many unnecessary variables
    # are created. One way to reduce this is to put all containers in a
    # canonical form (shortest side, middle side, longest side). Nevertheless,
    # if the variables are created but not used in constraints, I think they
    # do not cause major performance problems (to be verified).
    x[1:m, 1:6, 1:C, 1:maximum(X), 1:maximum(Y), 1:maximum(Z)], Bin
    e[1:C], Bin
  end

  @objective(model, Min, sum(e[j] * v[j] for j = 1:C))

  # Constraint (13) -- no overlap.
  for j = 1:C, s = 1:X[j], t = 1:Y[j], u = 1:Z[j]
    @constraint(model, sum(
      x[i,g,j,p,q,r] for i = 1:m, g = 1:6, p = max(1,s-l[i,g]+1):s, 
        q = max(1,t-w[i,g]+1):t, r = max(1,u-h[i,g]+1):u
    ) <= e[j])
  end

  # Constraint (14) -- every box needs to be packed (and it needs to be in a
  # valid non-protuding-outside-the-container position).
  for i = 1:m
    @constraint(model, sum(
      x[i,g,j,p,q,r] for g = 1:6, j = 1:C, p = 1:(X[j]-l[i,g]+1),
        q = 1:(Y[j]-w[i,g]+1), r = 1:(Z[j]-h[i,g]+1)
    ) == 1)
  end

  model
end

"""
    extract_attributions(model) :: Dict{Int64, Vector{Int64}}

Returns a dictionary with the index of the containers as keys and the index
of the boxes inside that container as values.
"""
function extract_attributions(model) :: Dict{Int64, Vector{Int64}}
  # TODO: this is slower than it could be. If the dimensions of every
  # box and each container were passed to this function, then the loops
  # could iterate a lot less unused variables. Nevertheless, 99% is 
  # spent on solving the model, retrieving the solution is not 
  # performance critical. More than that, is possible that getvalue
  # already copy the Matrix, so the damage is already done.
  x = getvalue(model[:x])
  m = size(x, 1) # number of boxes
  # size(x, 2) will be always six
  C = size(x, 3) # number of containers
  X = size(x, 4)
  Y = size(x, 5)
  Z = size(x, 6)

  d = Dict{Int64, Vector{Int64}}()
  # builds the dictionary from container index to box indexes
  for i = 1:m
    for g = 1:6, j = 1:C, p = 1:X, q = 1:Y, r = 1:Z
      if isone(round(x[i,g,j,p,q,r]))
        if haskey(d, j)
          push!(d[j], i)
        else
          d[j] = [i]
        end
        break # we already found where the box is so we can go to the next box
      end
    end
  end
  d
end

"""
    extract_positions(model) :: NTuple{3, Vector{Int64}}

Returns three arrays (x, y, z) of coordinates where the boxes were placed.
"""
function extract_positions(model) :: NTuple{3, Vector{Int64}}
  # TODO: see extract_attributions comment about the performance.
  x = getvalue(model[:x])
  m = size(x, 1) # number of boxes
  # size(x, 2) will be always six
  C = size(x, 3) # number of containers
  X = size(x, 4)
  Y = size(x, 5)
  Z = size(x, 6)

  s = Vector{Int64}(undef, m)
  t = Vector{Int64}(undef, m)
  u = Vector{Int64}(undef, m)
  # find and save the position of every box
  for i = 1:m
    for g = 1:6, j = 1:C, p = 1:X, q = 1:Y, r = 1:Z
      if isone(round(x[i,g,j,p,q,r]))
        # we change the coordinate system from base one to base zero here
        s[i], t[i], u[i] = p-1, q-1, r-1
        break # we already found where the box is so we can go to the next box
      end
    end
  end
  
  (s, t, u)
end

function orientate(g, l, w, h)
  if     g == 1; (l, w, h)
  elseif g == 2; (l, h, w)
  elseif g == 3; (w, l, h)
  elseif g == 4; (w, h, l)
  elseif g == 5; (h, l, w)
  elseif g == 6; (h, w, l)
  else
    error("The g parameter of Kurpel2018.orientate should stay between " *
          "1 and 6 (both inclusive) but it was " * string(g) * ".")
  end
end

"""
    extract_oriented_sizes(model, pqr) :: NTuple{3, Vector{Int64}}

Returns the oriented sizes of the boxes in the solved model.
"""
function extract_oriented_sizes(model, LWH) :: NTuple{3, Vector{Int64}}
  # TODO: see extract_attributions comment about the performance.
  x = getvalue(model[:x])
  m = size(x, 1) # number of boxes
  # size(x, 2) will be always six
  C = size(x, 3) # number of containers
  X = size(x, 4)
  Y = size(x, 5)
  Z = size(x, 6)

  L, W, H = LWH

  l = Vector{Int64}(undef, m)
  w = Vector{Int64}(undef, m)
  h = Vector{Int64}(undef, m)
  # discover which orientation every box is and give the oriented dimensions
  for i = 1:m
    for g = 1:6, j = 1:C, p = 1:X, q = 1:Y, r = 1:Z
      if isone(round(x[i,g,j,p,q,r]))
        l[i], w[i], h[i] = orientate(g, L[i], W[i], H[i])
        break # we already found where the box is so we can go to the next box
      end
    end
  end

  (l, w, h)
end

"""
    extract_solution(model, LWH)

Returns the triple of the results of extract\\_{attributions, positions,
oriented\\_sizes} over the solved model.

This is faster than executing the three functions separately.
"""
function extract_solution(model, LWH)
  # TODO: see extract_attributions comment about the performance.
  x = getvalue(model[:x])
  m = size(x, 1) # number of boxes
  # size(x, 2) will be always six
  C = size(x, 3) # number of containers
  X = size(x, 4)
  Y = size(x, 5)
  Z = size(x, 6)

  L, W, H = LWH

  # for some container, which boxes are inside it
  d = Dict{Int64, Vector{Int64}}()
  # coordinate of the box inside the container
  s = Vector{Int64}(undef, m)
  t = Vector{Int64}(undef, m)
  u = Vector{Int64}(undef, m)
  # the dimensions of the items when already rotated
  l = Vector{Int64}(undef, m)
  w = Vector{Int64}(undef, m)
  h = Vector{Int64}(undef, m)

  for i = 1:m
    for g = 1:6, j = 1:C, p = 1:X, q = 1:Y, r = 1:Z
      if isone(round(x[i,g,j,p,q,r]))
        if haskey(d, j)
          push!(d[j], i)
        else
          d[j] = [i]
        end
        # we change the coordinate system from base one to base zero here
        s[i], t[i], u[i] = p-1, q-1, r-1
        l[i], w[i], h[i] = orientate(g, L[i], W[i], H[i])
        break # we already found where the box is so we can go to the next box
      end
    end
  end

  (d, (s, t, u), (l, w, h))
end

end # module Kurpel2018

