#include <cstdint> // for uint*_t
#include <algorithm> // for sort
#include <vector>

using namespace std;

#include "exact_corner_points.hpp"
using namespace hbd;

#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#include "catch.hpp"

// I did not bother to discover how to make templates work with catch, then
// this is a dirty workaround while the test cases are not properly
// templated.
typedef uint16_t S;
typedef uint16_t BI;
typedef uint16_t CPI;

template <typename S>
vector<array<S, 6>> scratch_sextuple_set(
  triple<vector<S>> &b_coordinates,
  triple<vector<S>> &b_dimensions,
  triple<S> container_dims
) {
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

  return scratch_free_cps;
}

template <typename S, typename BI, typename CPI>
vector<array<S, 6>> delta_sextuple_set(
  const triple<vector<S>> &b_coordinates,
  const triple<vector<S>> &b_dimensions,
  triple<S> container_dims
) {
  space<S, BI, CPI> s(container_dims);
  for (BI i = 0; i < b_coordinates[0].size(); ++i) {
    CPI cp_ix = 0;
    for (; cp_ix < s.cp[0].size(); ++cp_ix) {
      if (
        s.cp[0][cp_ix] == b_coordinates[0][i] &&
        s.cp[1][cp_ix] == b_coordinates[1][i] &&
        s.cp[2][cp_ix] == b_coordinates[2][i]
      ) break;
    }

    // The corner point given in the test data must exist.
    assert(cp_ix < s.cp[0].size());

    const triple<S> b_dim{
      b_dimensions[0][i], b_dimensions[1][i], b_dimensions[2][i]
    };
    s.place(b_dim, cp_ix);
  }
  
  vector<array<S, 6>> delta_free_cps;
  for (CPI cp_ix = 0; cp_ix < s.cp[0].size(); ++cp_ix) {
    if (!s.clog_nb[cp_ix]) {
      delta_free_cps.push_back(array<S, 6>{
        s.cp[0][cp_ix],
        s.cp[1][cp_ix],
        s.cp[2][cp_ix],
        s.cp_min_b_after[0][cp_ix] - s.cp[0][cp_ix],
        s.cp_min_b_after[1][cp_ix] - s.cp[1][cp_ix],
        s.cp_min_b_after[2][cp_ix] - s.cp[2][cp_ix]
      });
    }
  }
  sort(delta_free_cps.begin(), delta_free_cps.end());

  return delta_free_cps;
}

void are_scratch_and_delta_same_as(
  triple<vector<S>> &b_coordinates,
  triple<vector<S>> &b_dimensions,
  triple<S> container_dims,
  vector<array<S, 6>> &expected_cps
) {
  sort(expected_cps.begin(), expected_cps.end());
  WHEN("comparing the values computed by scratch and delta") {
    const vector<array<S, 6>> scratch_free_cps = scratch_sextuple_set<S>(
      b_coordinates, b_dimensions, container_dims
    );
    const vector<array<S, 6>> delta_free_cps = delta_sextuple_set<S, BI, CPI>(
      b_coordinates, b_dimensions, container_dims
    );

    THEN("the scratch method return the expected free cps") {
      REQUIRE(scratch_free_cps == expected_cps);
    }
    THEN("the delta method return the expected free cps") {
      REQUIRE(delta_free_cps == expected_cps);
    }
  }
}

// I think I probably should learn how to write good scenarios.
SCENARIO("computing free cps from scratch or delta works") {
  // Will use the same container dimensions for all tests to simplify.
  triple<S> container_dims{100, 200, 400};
  GIVEN("a single-box instance") {
    triple<vector<S>> b_coordinates{{{1}, {1}, {1}}};
    triple<vector<S>> b_dimensions{{{50}, {100}, {200}}};
    vector<array<S, 6>> expected_cps{
      {51,   1,   1, 1, 1, 1},
      { 1, 101,   1, 1, 1, 1},
      { 1,   1, 201, 1, 1, 1}
    };
    are_scratch_and_delta_same_as(
      b_coordinates, b_dimensions, container_dims, expected_cps
    );
  }
  GIVEN("a single box using all space on all dimensions") {
    triple<vector<S>> b_coordinates{{{1}, {1}, {1}}};
    triple<vector<S>> b_dimensions{{{100}, {200}, {400}}};
    vector<array<S, 6>> expected_cps{};
    are_scratch_and_delta_same_as(
      b_coordinates, b_dimensions, container_dims, expected_cps
    );
  }
  GIVEN("a single box using all space on a single dimension") {
    triple<vector<S>> b_coordinates{{{1}, {1}, {1}}};
    triple<vector<S>> b_dimensions{{{100}, {100}, {200}}};
    vector<array<S, 6>> expected_cps{
      {1, 101,   1, 1, 1, 1},
      {1,   1, 201, 1, 1, 1}
    };
    are_scratch_and_delta_same_as(
      b_coordinates, b_dimensions, container_dims, expected_cps
    );
  }
  GIVEN("four boxes (all first three cps used) nine cps") {
    triple<vector<S>> b_coordinates{{
      {1, 51, 1, 1}, {1, 1, 101, 1}, {1, 1, 1, 151}
    }};
    triple<vector<S>> b_dimensions{{
      {50, 20, 30, 40}, {100, 40, 50, 60}, {150, 70, 80, 90}
    }};
    vector<array<S, 6>> expected_cps{
       { 1,  61, 151, 1, 1, 1},
       {41,   1, 151, 1, 1, 1},
       {51,   1,  71, 1, 1, 1},
       { 1, 101,  81, 1, 1, 1},
       { 1,   1, 241, 1, 1, 1},
       {31, 101,   1, 1, 1, 1},
       {51,  41,   1, 1, 1, 1},
       { 1, 151,   1, 1, 1, 1},
       {71,   1,   1, 1, 1, 1}
       /*{ 1,  70, 160,  2,  71, 161},
       {50,   1, 160, 51,   2, 161},
       {60,   1,  80, 61,   2,  81},
       { 1, 110,  90,  2, 111,  91},
       { 1,   1, 250,  2,   2, 251},
       {40, 110,   1, 41, 111,   2},
       {60,  50,   1, 61,  51,   2},
       { 1, 160,   1,  2, 161,   2},
       {80,   1,   1, 81,   2,   2}*/
    };
    are_scratch_and_delta_same_as(
      b_coordinates, b_dimensions, container_dims, expected_cps
    );
  }
  GIVEN("tricky, 2 boxes, 3 cps, 1 concrete, 1 not, and 1 both") {
    triple<vector<S>> b_coordinates{{{1, 51}, {1, 1}, {1, 1}}};
    triple<vector<S>> b_dimensions{{{50, 50}, {100, 100}, {200, 300}}};
    vector<array<S, 6>> expected_cps{
      {1,   1, 201,  1, 1, 1},
      {1,   1, 301, 51, 1, 1},
      {1, 101,   1,  1, 1, 1},
      {1, 101,   1, 51, 1, 1}
    };
    are_scratch_and_delta_same_as(
      b_coordinates, b_dimensions, container_dims, expected_cps
    );
  }
  GIVEN("") {
    triple<vector<S>> b_coordinates{{
      {1, 11, 11, 51,   1,   1},
      {1,  1,  1,  1, 101,   1},
      {1,  1, 11,  1,   1, 151}
    }};
    triple<vector<S>> b_dimensions{{
      { 10,  40,  40, 10, 30, 20},
      {100, 100,  10, 50, 10, 60},
      {150,  10, 140, 70, 60, 10}
    }};
    vector<array<S, 6>> expected_cps{
      { 1,   1, 161,  1,  1,   1},
      { 1,  61, 151,  1,  1,   1},
      { 1, 101,  61,  1,  1,   1},
      { 1, 111,   1,  1,  1,   1},
      {11,  11,  11,  1,  1,   1}, 
      {11,  11,  61,  1, 91,   1},
      {11,  11,  71, 41,  1,   1}, 
      {11,  51,  11, 41,  1,   1}, 
      {11,  51,  61, 41, 51,   1}, 
      {11,  61,  11,  1,  1, 141},
      {11,  61,  61,  1, 41,  91},
      {21,   1, 151,  1,  1,   1},
      {21,  11,  11,  1,  1, 141},
      {21,  11,  61,  1, 91,  91},
      {21,  11,  71, 31,  1,  81}, 
      {21,  51,  11, 31,  1, 141},
      {21,  51,  61, 31, 51,  91},
      {31,  11,  11,  1, 91,   1},
      {31,  51,  11, 21, 51,   1}, 
      {31, 101,   1,  1,  1,   1},
      {51,   1,  71,  1,  1,   1}, 
      {51,  51,   1,  1,  1,   1}, 
      {61,   1,   1,  1,  1,   1} 
    };
    are_scratch_and_delta_same_as(
      b_coordinates, b_dimensions, container_dims, expected_cps
    );
  }
}

