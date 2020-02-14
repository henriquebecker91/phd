module ExtremumPlanes

const SizeType = UInt16
const VolumeType = UInt32

abstract type RectangularCuboid end

function volume(a :: RectangularCuboid) :: VolumeType
  VolumeType(a.l) * VolumeType(a.w) * VolumeType(a.h)
end

# let us start with no rotation to keep things simple
struct Box <: RectangularCuboid
  l :: SizeType
  w :: SizeType
  h :: SizeType

  function Box(l, w, h)
    @assert l > 0
    @assert w > 0
    @assert h > 0

    new(l, w, h)
  end
end

struct Point
  x :: SizeType
  y :: SizeType
  z :: SizeType

  Point(x, y, z) = new(x, y, z)
end

struct Container <: RectangularCuboid
  l :: SizeType
  w :: SizeType
  h :: SizeType

  boxes :: Vector{Box}
  corners :: Vector{Point}

  Container(l, w, h, boxes = [], corners = [Point(0, 0, 0)]) = new(l, w, h, boxes, corners)
end

end # module ExtremumPlanes

