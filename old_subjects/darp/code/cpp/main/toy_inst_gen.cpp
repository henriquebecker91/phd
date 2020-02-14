#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <cassert>
#include "../lib/darp_instances.hpp"

using namespace std;

template <typename P, typename V>
struct custom_int_distribution {
  std::vector<P> ps;
  std::vector<V> vs;
  std::uniform_int_distribution<P> d;
  
  custom_int_distribution(
    const std::vector<P> &ps_,
    const std::vector<V> &vs_
  ) : ps(ps_), vs(vs_) {
    assert(ps.size() == vs.size());
    for (P i = 1; i < ps.size(); ++i) ps[i] += ps[i - 1];
    d = uniform_int_distribution<P>(1, ps[ps.size() - 1]);
    ps.pop_back();
  }

  template<typename URNG>
  V operator()(URNG &g) {
    const P p = d(g);
    for (P i = 0; i < ps.size(); ++i) if (p <= ps[i]) return vs[i];

    return vs.back();
  }
};

int main (void)
{
  mt19937_64 rng(0);
  const vector<uint8_t> ps1{95, 5};
  const vector<uint8_t> vs1{4, 7};
  custom_int_distribution<uint8_t, uint8_t> dvc(ps1, vs1);
  const vector<uint8_t> ps2{60, 30, 8, 2};
  const vector<uint8_t> vs2{1, 2, 3, 4};
  custom_int_distribution<uint8_t, uint8_t> drp(ps2, vs2);
  /*vector<uint16_t> qt{0, 0, 0, 0, 0};
  for (uint16_t i = 0; i < 1000; ++i) {
    ++qt[cid(rng)];
  }
  for (const uint16_t &v : qt) {
    cout << v << endl;
  }*/

  coord_darp_inst_t<
    uint16_t,
    uint16_t,
    uint64_t,
    uint8_t> inst(
    3,
    50,
    1000,
    1000,
    rng,
    dvc,
    drp
  );
  //cout << inst.to_s();
  cout << inst.to_dot();

	return EXIT_SUCCESS;
}

