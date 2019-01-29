#include "extremum_planes.hpp"
using namespace hbm;

#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#include "catch.hpp"

#include <stdint.h> // for uint32_t

typedef uint16_t length_t;
typedef uint32_t volume_t;

SCENARIO("containers have CP's where boxes can be placed") {
  GIVEN("an empty container") {
    
    container_t<length_t> c(10, 20, 30);

    REQUIRE(c.cps.front() == cp_t<length_t>(0, 0, 0));
    REQUIRE(c.cps.cbegin() == --c.cps.cend());

    WHEN("") {
      THEN("") {
        REQUIRE(1 == 1);
      }
    }
  }
}

