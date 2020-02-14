using CornerPoints

function transpack(xs)
  n = (unique ∘ map)(length, xs)
  @assert (isone ∘ length)(n) # all vectors inside xs have the same length
  n = n[1] # convert from size one array to scalar
  ys = []
  for i = 1:n
    push!(ys, map(a->a[i], xs))
  end
  ys
end

# st: set of tuples
# ta: tuple of arrays
# given Set{(1, 2), (3, 4)} and ([1, 3], [2, 4]) it will give true, i.e., it
# checks if the set of tuples and the tuple of arrays store the same data,
# disregarding the order between the elements themselves (the set will be
# always ordered and the ta can be any permutation of the same elements).
function same(dt, ta)
  stta = Set(transpack(ta))
  dt == stta
end

@testset "CornerPoints" begin
  @testset "compute_free_cps" begin
    LWH = (10, 20, 30)
    @testset "single box smaller than container in all dimensions" begin
      xyzpqr = ([1], [1], [1], [5], [10], [15])
      answer = Set([
        (1,  1, 16, 1, 1, 1, 2, 3, 1),
        (1, 11,  1, 1, 1, 1, 2, 1, 4),
        (6,  1,  1, 1, 1, 1, 1, 3, 4)
      ])
      @test same(answer, compute_free_cps(xyzpqr, LWH))
    end
    @testset "single box using all space in one dimension" begin
      xyzpqr = ([1], [1], [1], [10], [10], [15])
      answer = Set([
        (1,  1, 16, 1, 1, 1, 2, 3, 1),
        (1, 11,  1, 1, 1, 1, 2, 1, 4)
      ])
      @test same(answer, compute_free_cps(xyzpqr, LWH))
    end
    @testset "single box using all space in all dimensions" begin
      xyzpqr = ([1], [1], [1], [10], [20], [30])
      answer = Set()
      @test same(answer, compute_free_cps(xyzpqr, LWH))
    end
    @testset "single box using all space in all dimensions" begin
      xyzpqr = ([1], [1], [1], [10], [20], [30])
      answer = Set()
      @test same(answer, compute_free_cps(xyzpqr, LWH))
    end
    @testset "four boxes (all first three cps used) nine cps" begin
      LWH = (10, 20, 30)
      xyzpqr = (
        [ 1, 6,  1,  1],
        [ 1, 1, 11,  1],
        [ 1, 1,  1, 16],
        [ 5, 2,  3,  4],
        [10, 4,  5,  6],
        [15, 7,  8,  9]
      )
      answer = Set([
       (1,  7, 16, 1, 1, 1, 5, 4, 1),
       (5,  1, 16, 1, 1, 1, 4, 6, 1),
       (6,  1,  8, 1, 1, 1, 1, 6, 2),
       (1, 11,  9, 1, 1, 1, 5, 1, 3),
       (1,  1, 25, 1, 1, 1, 5, 6, 4),
       (4, 11,  1, 1, 1, 1, 3, 1, 7),
       (6,  5,  1, 1, 1, 1, 1, 2, 7),
       (1, 16,  1, 1, 1, 1, 5, 3, 7),
       (8,  1,  1, 1, 1, 1, 2, 6, 7)
      ])
      @test same(answer, compute_free_cps(xyzpqr, LWH))
    end
    @testset "tricky, 2 boxes, 3 cps, 1 concrete, 1 not, and 1 both" begin
      xyzpqr = ([1, 6], [1, 1], [1, 1], [5, 5], [10, 10], [15, 20])
      # The first box is half each dimension. The second box fill up the first
      # dimension, leaves the same gap as the first box in the second
      # dimension, and exceeds the first box in the third dimension. This way,
      # the second box display three different behaviours: it has no corner
      # point with height and width one (because lenght was all filled up); the
      # corner point of lenght and height one is the same as the one of the
      # first box but with a min-box-length of six; the corner point with
      # length and width one is above the one of the first box, but has
      # min-box-length of six.
      answer = Set([
        (1, 1, 16, 1, 1, 1, 3, 4, 1)
        (1, 1, 21, 6, 1, 1, 3, 4, 2)
        (1, 11, 1, 1, 1, 1, 3, 1, 5)
        (1, 11, 1, 6, 1, 1, 3, 2, 5)
      ])

      @test same(answer, compute_free_cps(xyzpqr, LWH))
    end
    @testset "the 'entirely in the air' corner point" begin
      xyzpqr = ([1, 2, 2, 6, 1, 1], [1, 1, 1, 1, 11, 1], [1, 1, 2, 1, 1, 16], [1, 4, 4, 1, 3, 2], [10, 10, 1, 5, 1, 6], [15, 1, 14, 7, 6, 1])
      answer = Set([
        (1, 1, 17, 1, 1, 1, 7, 8, 6),
        (1, 7, 16, 1, 1, 1, 7, 6, 1),
        (1, 11, 7, 1, 1, 1, 7, 1, 5),
        (1, 12, 1, 1, 1, 1, 7, 5, 9),
        (2, 2, 2, 1, 1, 1, 1, 3, 2),
        (2, 2, 7, 1, 10, 1, 1, 3, 5),
        (2, 2, 8, 5, 1, 1, 1, 3, 4),
        (2, 6, 2, 5, 1, 1, 1, 4, 2),
        (2, 6, 7, 5, 6, 1, 1, 4, 5),
        (2, 7, 2, 1, 1, 15, 1, 6, 2),
        (2, 7, 7, 1, 5, 10, 1, 6, 5),
        (3, 1, 16, 1, 1, 1, 6, 8, 3),
        (3, 2, 2, 1, 1, 15, 6, 3, 2),
        (3, 2, 7, 1, 10, 10, 6, 3, 5),
        (3, 2, 8, 4, 1, 9, 6, 3, 4),
        (3, 6, 2, 4, 1, 15, 6, 4, 2),
        (3, 6, 7, 4, 6, 10, 6, 4, 5),
        (4, 2, 2, 1, 10, 1, 5, 3, 2),
        (4, 6, 2, 3, 6, 1, 5, 4, 2),
        (4, 11, 1, 1, 1, 1, 5, 2, 9),
        (6, 1, 8, 1, 1, 1, 3, 8, 4),
        (6, 6, 1, 1, 1, 1, 2, 4, 9),
        (7, 1, 1, 1, 1, 1, 4, 8, 9)
      ])

      @test same(answer, compute_free_cps(xyzpqr, LWH))
    end
  end
end

