#include <cstdlib>

#include "../lib/darp_instances.hpp"

template <
  template <typename, typename, typename, typename, typename> class O,
  typename V,
  typename L,
  typename P,
  typename C,
  typename T>
int templated_main(int argc, char** argv) {
  signal(SIGINT, &catch_function);
  O<V, L, P, C, T> inst(std::cin);
  inst.debug_output();

  C bkv;
  std::vector<set_map_t<L, L>*> bkv_sol;
  test_all_mv_schedules<O, V, L, P, C, T>(inst, bkv_sol, bkv);

  return EXIT_SUCCESS;
}

int main(int argc, char** argv) {
  return templated_main<CL2003_darp_inst_wt, short, short, short, double, double>(argc, argv);
}

