#include "extremum_planes.hpp"
using namespace hbm;

#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#include "catch.hpp"

#include <stdint.h> // for uint32_t

typedef uint16_t length_t;
typedef uint32_t volume_t;

SCENARIO("both container and box types have volume") {
  GIVEN("a box and a container with large dimensions") {
    box_t<length_t, volume_t> b(100, 200, 300);
    container_t<length_t, volume_t> c(300, 200, 100);

    WHEN("computing a volume too big for the L type") {
      volume_t bv = b.volume();
      volume_t cv = c.volume();

      THEN("we yet get the right values") {
        REQUIRE(bv == 6000000)
        REQUIRE(cv == 6000000)
      }
    }
  }
}

SCENARIO("a space has CPs where boxes can be placed") {
  GIVEN("an empty space") {
    
    space<length_t, volume_t> s(10, 20, 30);

    REQUIRE(c.open_cps.front() == cp<length_t>(0, 0, 0));
    REQUIRE(c.open_cps.cbegin() == --c.open_cps.cend());

    auto origin_cp = c.open_cps.begin();

    REQUIRE(origin_cp->free());

    WHEN("placing a box with a quarter of the volume") {
      box_t b(5, 10, 15);

      s.place(b, origin_cp);

      THEN("the first CP is now marked as used") {
        REQUIRE(origin_cp->used());
        REQUIRE(s.used_cps.begin() == origin_cp);
      }

      THEN("three new open CPs are created") {
        auto fcps_cbegin = s.free_cps.cbegin();
        auto fcps_cend = s.free_cps.cend();
        REQUIRE(count_if(fcps_cbegin, fcps_cend, [](auto &cp) {
          return true;
        }) == 3)
        REQUIRE(find_if(fcps_cbegin, fcps_cend, [](auto &cp) {
          return cp.x == 5 && cp.y == 0 && cp.z == 0 && cp.free();
        }) != fcps_cend);
        REQUIRE(find_if(fcps_cbegin, fcps_cend, [](auto &cp) {
          return cp.x == 0 && cp.y == 10 && cp.z == 0 && cp.free();
        }) != fcps_cend);
        REQUIRE(find_if(fcps_cbegin, fcps_cend, [](auto &cp) {
          return cp.x == 0 && cp.y == 0 && cp.z == 15 && cp.free();
        }) != fcps_cend);
      }

      THEN("") {

      }
    }

    WHEN("placing a box larger than the container in some dimension") {
      box_t b(5, 10, 31);

      s.place(b, origin_cp);
    }
  }

  /* place a box smaller than the container
   * place a box greater than the container in one dimension (times three)
   * place a box greater than all container dimensions
   * place a box smaller than all container dimensions in a cp where
   *   it ends up protruding from the container (times three)
   * place a box in a way it overlaps with one other box
   * place a box in a way it overlaps with two other boxes
   * place a box in a way it overlaps with some box and protrudes
   *   the container
   * place a box in a cp with a projected box plane, in a way that breaks
   *   the cannot-push-left-back-bottom assumption
   * place a box in a cp with a projected box plane, in a way that respects
   *   the cannot-push-left-back-bottom assumption
   * place a box breaking the three rules: protrudes the container,
   *   overlaps with a box, and can be pushed left, back, or bottom
   *
   * remove a box that can be removed: check box list and cp list
   * test remove box that cannot be removed: check dependent cps and
   *   dismissed cps
   *
   * place box and check if the right cps were generated
   *   one single box and the three usual cps
   *   if a box face touch the container wall the related cp does not exist
   *   put a second box of the same size as the first in each possible cp
   *     and create only one cp (the other two are the same)
   *
   */
}

