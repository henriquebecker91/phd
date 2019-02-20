#ifndef HBM_EXTREMUM_PLANES_HPP
#define HBM_EXTREMUM_PLANES_HPP

#include <stdint> // for uint8_t in enum struct
#include <ostream> // for declaring how structs should be printed
#include <list>
#include <vector>

namespace hbm {
  // A corner point is always the meeting of three planes. However, two
  // distinct boxes can have overlapping box projections, and if they
  // are overlapping in the corner point, then removing one of the boxes will
  // not remove the corner point (even if it is associated to a plane that
  // is associated to that box).

  // Notation rules: a class will have the "_t" suffix if their objects are
  // immutable data given as parameters for the problem. Other class will
  // contain a mutable representation of the object for manipulation by the
  // algorithm. Such mutable representation will contain a copy (or reference)
  // of the immutable box data given by the instance.

  // data structure guideline: the main data structure 'space' will make use
  // of double-linked lists (referred to as 'lists' from now on), and the code
  // will rely on having O(1) operations removing or adding nodes.
  // Consequently, it is a good thing to have medium/large structs that will
  // be inside such nodes, as this reduces the pointer overhead.

  template <typename L, typename V>
  struct rect_cuboid_t {
    L l; L w; L h;

    rect_cuboid_t(L l, L w, L h) : x(l), y(w), z(h) {}

    // for now volume is a procedure, lets check in the future if its the
    // best to compute it in the construction and store it
    inline V volume(void) {
      // multiplication has left to right associativity so this code is safe
      return static_cast<V>(l) * w * h;
    }
  };

  // To keep things simple, we start with no rotation nor any other box
  // characteristic (profit, load bearing, etc...). 
  template <typename L, typename V>
  struct box_t : public rect_cuboid_t<L, V> {
    box_t(L l, L w, L h) : rect_cuboid_t(l, w, h) {}

    friend std::ostream& operator<<(std::ostream& os, const box_t &b) {
      return os << "box_t{" << b.l << "," << b.w << "," << b.h << "}";
    }
  };

  template <typename L, typename V>
  struct container_t : public rect_cuboid_t<L, V> {
    container_t(L l, L w, L h) : rect_cuboid_t(l, w, h) {}

    friend std::ostream& operator<<(std::ostream& os, const box_t &b) {
      return os << "container_t{" << b.l << "," << b.w << "," << b.h << "}";
    }
  };

  /*enum struct plane_base : uint8_t {
    container_wall,
    box_face,
    box_projection
  };
  enum struct plane_axis : uint8_t {
    xplane,
    yplane,
    zplane
  };

  template <typename L>
  struct plane {
    plane_base base;
    plane_axis axis;
    std::vector<cp_ref> associated_cps;
  };*/

  // Forward declare cp, as placed_box will have a cp_ref and cp will have
  // a placed_box_ref.
  template <typename L, typename V>
  struct cp;
  typedef std::list<cp<L, V>>::iterator cp_ref;

  template <typename L, typename V>
  struct placed_box {
    cp_ref<L> cp;
    box_t<L, V> b;
  };
  typedef std::list<placed_box<L, V>>::iterator placed_box_ref;

  // corner point, L = lenght type
  // This class is mostly for internal use. The space structure will
  // create the cp objects, and return references to them. Such references
  // are passed to other methods of the space structure.
  template <typename L, typename V>
  struct cp {
    L x; L y; L z;
    bool ff; // free flag, both free and used are inline methods
    // does cp really need a placed_box_ref? in which situation the caller
    // will have a used cp and would want to know which box is placed on it?

    std::list<placed_box_ref> xplane;
    std::list<placed_box_ref> yplane;
    std::list<placed_box_ref> zplane;

    cp(L x, L y, L z) : x(x), y(y), z(z), free(true) {}

    inline bool free(void) { return  ff; }
    inline bool used(void) { return !ff; }

    inline bool operator==(const cp_t& o) const {
      return x == o.x && y == o.y && z == o.z;
    }

    friend std::ostream& operator<<(std::ostream& os, const cp_t& cp) {
      return os << "cp{" << cp.x << "," << cp.y << "," << cp.z << "}";
    }
  };

  template <typename L, typename V>
  struct space {
    container_t<L, V> ctype;
    std::list<cp<L, V>> free_cps; // corner points without a box placed on them
    std::list<cp<L, V>> used_cps; // corner points with a box placed on them
    std::list<placed_box<L, V>> boxes; // packed boxes, same size as used_cps

    space(container_t<L, V> ctype) : ctype(ctype) {}


    pair<placed_box_ref, std::vector<cp<L, V>>> place(
      const box_t<L> &b, const cp_ref &cp
    ) {
      
    }
    
  };

}

#endif // HBM_EXTREMUM_PLANES_HPP

