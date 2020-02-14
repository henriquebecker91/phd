#include "chen1995_model.hpp"
#include <cstdlib>

#ifndef HBD_PRINT_VAR
  #define HBD_PRINT_VAR(var) out << (#var ": ") << var << std::endl
#endif

#ifndef HBD_PRINT_POS
  #define HBD_PRINT_POS(a, ix) out << (#a "[") << ix << ("]: ") << \
  a[ix] << std::endl
#endif

#ifndef HBD_PRINT_POS2
  #define HBD_PRINT_POS2(a, ix, ix2) out << (#a "[") << ix << "][" << \
  ix2 << "]: " << a[ix][ix2] << std::endl
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
  bool *lx = new bool[N], *ly = new bool[N], *lz = new bool[N],
       *wx = new bool[N], *wy = new bool[N], *wz = new bool[N],
       *hx = new bool[N], *hy = new bool[N], *hz = new bool[N];
  bool **a = new bool*[N], **b = new bool*[N], **c = new bool*[N],
       **d = new bool*[N], **e = new bool*[N], **f = new bool*[N];

  for (Q i = 0; i < N - 1; ++i) {
    a[i] = new bool[N]; b[i] = new bool[N]; c[i] = new bool[N]; 
    d[i] = new bool[N]; e[i] = new bool[N]; f[i] = new bool[N]; 
  }

  hbm::chen1995<Q, C>(
    N, m, s, n, p, q, r, &L, &W, &H, x, y, z, lx, ly, lz, wx, wy, wz,
    hx, hy, hz, a, b, c, d, e, f
  );
  hbm::chen1995_check<Q, C>(
    N, m, s, p, q, r, &L, &W, &H, x, y, z, lx, ly, lz, wx, wy, wz, hx, hy, hz
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

  // freeing up memory in reverse order
  for (Q i = 0; i < N - 1; ++i) {
    delete[] a[i]; delete[] b[i]; delete[] c[i]; 
    delete[] d[i]; delete[] e[i]; delete[] f[i]; 
  }
  delete[] a; delete[] b; delete[] c; delete[] d; delete[] e; delete[] f;
  delete[] lx; delete[] ly; delete[] lz; 
  delete[] wx; delete[] wy; delete[] wz; 
  delete[] hx; delete[] hy; delete[] hz; 
  delete[] x; delete[] y; delete[] z; 
  delete[] n;
  for (Q i = 0; i < N; ++i) delete[] s[i];
  delete[] s;

  return EXIT_SUCCESS;
}


// Solves an UKP instance withe the gurobi solver. Probably print lots
// of extra information in addition to the optimal solution.
// Takes the name of a file in the ".ukp" format. Other options
// should be consulted at gurobi_ukp_model.hpp.
int main(int argc, char** argv) {
  auto &out = std::cout;
  HBD_PRINT_VAR(HBM_GIT_HEAD_AT_COMPILATION);
  return templated_main<short, short>(argc, argv);
}

