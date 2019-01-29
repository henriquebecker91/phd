#ifndef HBM_EXTREMUM_PLANES_HPP
#define HBM_EXTREMUM_PLANES_HPP

#include <list>

namespace hbm {
  // A corner point is always the meeting of three planes. However, two
  // distinct boxes can have overlapping box projections, and if they
  // are overlapping in the corner point, then removing one of the boxes will
  // not remove the corner point (even if it is associated to a plane that
  // is associated to that box).

  // corner point type: L = lenght type
  template <typename L>
  struct cp_t {
    L x;
    L y;
    L z;

    cp_t(L x, L y, L z) : x(x), y(y), z(z) {}

    inline bool operator==(const cp_t &o) const {
      return x == o.x && y == o.y && z && o.z;
    }
  };

  template <typename L>
  struct container_t {
    L l;
    L w;
    L h;

    std::list<cp_t<L>> cps; // corner points

    container_t(L l, L w, L h) : l(l), w(w), h(h), cps(1, cp_t<L>(0, 0, 0)) { }
  };
}

#endif // HBM_EXTREMUM_PLANES_HPP

