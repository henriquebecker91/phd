extern "C" {
  #include "exact_corner_points_API_C.h"
}

#include "exact_corner_points.hpp"

typedef hbd::space<itype, ntype, ctype> ispace; // internal space type
typedef hbd::triple<itype> xyz;

#define DREF(s) reinterpret_cast<ispace*>(s)

extern "C" {
  space* new_space(itype l, itype w, itype h) {
    return reinterpret_cast<space*>(new ispace(xyz{l, w, h}));
  }
  void del_space(space *s) { delete DREF(s); }
  void pop(space *s) { DREF(s)->pop(); }

  void place(space *s, itype l, itype w, itype h, ctype cp_ix) {
      DREF(s)->place(xyz{l, w, h}, cp_ix);
  }

  boolean can_be_placed(space *s, itype l, itype w, itype h, ctype cp_ix) {
    return static_cast<boolean>(DREF(s)->can_be_placed(xyz{l, w, h}, cp_ix));
  }

  const ctype cpnb(space *s) {
    return static_cast<ctype>(DREF(s)->cp[0].size());
  }
  const itype* cpx(space *s) { return DREF(s)->cp[0].data(); }
  const itype* cpy(space *s) { return DREF(s)->cp[1].data(); }
  const itype* cpz(space *s) { return DREF(s)->cp[2].data(); }
  /*const ntype* cp_clog(space *s) { return DREF(s)->clog_nb.data(); }*/
}

