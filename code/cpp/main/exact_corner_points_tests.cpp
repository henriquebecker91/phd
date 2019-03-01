#include <cstdint> // for uint*_t
#include <algorithm> // for sort
#include <vector>

using namespace std;

#include "exact_corner_points.hpp"
using namespace hbm;

#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#include "catch.hpp"

// I did not bother to discover how to make templates work with catch, then
// this is a dirty workaround while I the test cases are not properly
// templated.
typedef uint16_t S;
typedef uint16_t BI;
typedef uint16_t CPI;

SCENARIO("corner points are the same if computed by scratch or delta") {
  GIVEN("a simple instance") {
    triple<vector<S>> b_coordinates{{1}, {1}, {1}};
    triple<vector<S>> b_dimensions{{5}, {10}, {15}}
    triple<S> container_dims{10, 20, 30};
    WHEN("comparing the values computed by scratch and delta") {
      array<triple<std::vector<S>>, 2> scratch_return = free_cps_from_scratch(
        b_coordinates, b_dimensions, container_dims
      );
      vector<array<S, 6>> scratch_free_cps;
      for (CPI cp_ix = 0; cp_ix < scratch_return[0][0].size(); ++cp_ix) {
        scratch_free_cps.push_back(array<S, 6>{
          scratch_return[0][0][cp_ix],
          scratch_return[0][1][cp_ix],
          scratch_return[0][2][cp_ix],
          scratch_return[1][0][cp_ix],
          scratch_return[1][1][cp_ix],
          scratch_return[1][2][cp_ix]
        });
      }
      sort(scratch_free_cps.begin(), scratch_free_cps.end());

      space<S> s(container_dims);
      for (BI i = 0; i < b_coordinates[0].size(); ++i) {
        CPI cp_ix = 0;
        for (; cp_ix < s.cp[0].size(); ++cp_ix) {
          if (
            s.cp[0][cp_ix] == b_coordinates[0][cp_ix] &&
            s.cp[1][cp_ix] == b_coordinates[1][cp_ix] &&
            s.cp[2][cp_ix] == b_coordinates[2][cp_ix]
          ) break;
        }

        array<triple<S>, 2> 
        // The corner point given in the test data must exist.
        assert(cp_ix < s.cp[0].size());
        const triple<S> b_dim {
          b_dimensions[0][i],
          b_dimensions[1][i],
          b_dimensions[2][i]
        };
        s.place(b_dim, cp_ix);
      }
      
      for (CPI cp_ix = 0; cp_ix < s.cp[0].size(); ++cp_ix) {
        
      }
      THEN("the two have the same set of free cps") {
        REQUIRE(scratch_free_cps == delta_free_cps)
      }
    }
  }
}

