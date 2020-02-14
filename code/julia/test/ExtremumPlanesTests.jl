using ExtremumPlanes

@testset "Basics (a.k.a. placing a single box)" begin
  c = Container(10, 10, 10)
  @test c.corners == [Point(0, 0, 0)]
  b1 = Box(5, 5, 5)
  place(b, c, 1)
end


