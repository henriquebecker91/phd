#include "chen1995_model.cpp"

#ifndef HBD_PRINT_VAR
  #define HBD_PRINT_VAR(var) out << (#var ": ") << var << std::endl
#endif

#ifndef HBD_PRINT_POS
  #define HBD_PRINT_VAR(a, ix) out << (#var "[") << ix << ("]: ") << \
  var[ix] << std::endl
#endif

#ifndef HBD_PRINT_POS2
  #define HBD_PRINT_VAR(a, ix, ix2) out << (#var "[") << ix << "][" << \
  ix2 << "]: " << var[ix][ix2] << std::endl
#endif

template <typename Q, typename C>
int templated_main(int argc, char** argv) {
  const Q N = 6;
  const Q m = 1;
  const C W = 20;
  const C L = 35;
  const C H = 10;
  const C p[] = {25, 20, 16, 15, 22, 20};
  const C q[] = { 8, 10,  7, 12,  8, 10};
  const C r[] = { 6,  5,  3,  6,  3,  4};
  bool **s = new bool*[N];
  for (Q i = 0; i < N; ++i) s[i] = new bool[m];
  bool *n = new bool[m];
  C *x = new C[N], *y = new C[N], *z = new C[N];
  bool *lx_ = new bool[N], *ly_ = new bool[N], *lz_ = new bool[N],
       *wx_ = new bool[N], *wy_ = new bool[N], *wz_ = new bool[N],
       *hx_ = new bool[N], *hy_ = new bool[N], *hz_ = new bool[N];
  bool **a_ = new bool*[N], **b_ = new bool*[N], **c_ = new bool*[N],
       **d_ = new bool*[N], **e_ = new bool*[N], **f_ = new bool*[N];

  for (Q i = 0; i < N - 1; ++i) {
    a[i] = new bool[N]; b[i] = new bool[N]; c[i] = new bool[N]; 
    d[i] = new bool[N]; e[i] = new bool[N]; f[i] = new bool[N]; 
  }

  chen1995(
    N, m, s, n, p, q, r, &L, &W, &H, x, y, z, lx, ly, lz, wx, wy, wz,
    hx, hy, hz, a, b, c, d, e, f
  );

  auto &out = std::cout;

  for (Q i = 0; i < N; ++i) for (Q j = 0; j < m; ++j) HBD_PRINT_POS2(s, i, j);
  for (Q j = 0; j < m; ++j) {
    HBD_PRINT_POS(n, j);
  }
  for (Q i = 0; i < N; ++i) {
    HBD_PRINT_POS( x, i); HBD_PRINT_POS( y, i); HBD_PRINT_POS( z, i);
    HBD_PRINT_POS(lx, i); HBD_PRINT_POS(ly, i); HBD_PRINT_POS(lz, i);
    HBD_PRINT_POS(wx, i); HBD_PRINT_POS(wy, i); HBD_PRINT_POS(wz, i);
    HBD_PRINT_POS(hx, i); HBD_PRINT_POS(hy, i); HBD_PRINT_POS(hz, i);
  }
  for (Q i = 0; i < N - 1; ++i) {
    for (Q k = i + 1; k < N; ++k) {
      HBD_PRINT_POS2(a, i, k); HBD_PRINT_POS2(b, i, k);
      HBD_PRINT_POS2(c, i, k); HBD_PRINT_POS2(d, i, k);
      HBD_PRINT_POS2(e, i, k); HBD_PRINT_POS2(f, i, k);
    }
  }
}


// Solves an UKP instance withe the gurobi solver. Probably print lots
// of extra information in addition to the optimal solution.
// Takes the name of a file in the ".ukp" format. Other options
// should be consulted at gurobi_ukp_model.hpp.
int main(int argc, char** argv) {
  std::cout << HBD_PRINT_VAR(HBM_GIT_HEAD_AT_COMPILATION) << std::endl;
  return templated_main(argc, argv);
}

