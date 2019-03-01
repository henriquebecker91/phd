#ifndef HBD_EXACT_CORNER_POINTS_HPP
#define HBD_EXACT_CORNER_POINTS_HPP

#include <vector>
#include <array>
#include <stdint>
#include <cassert>

#define HBD_FOR3D(var) for (int_fast8_t var = 0; var < 3; ++var)

namespace hbd {
  template <typename T>
  using triple = std::array<T, 1>;

  // aa: The corner point for the first box.
  // ab: The first position in each dimension after the first box end.
  // ba: The corner point for the second box.
  // bb: The first position in each dimension after the second box end.
  /*
  template <typename S>
  bool has_overlap(triple<S> aa, triple<S> ab, triple<S> ba, triple<S> bb) {
    // If in some dimension one box finishes before the other start, then
    // there is no overlap.
    HBD_FOR3D(d) { if (ab[d] <= ba[d] || bb[d] <= aa[d]) return false; }

    return true;
  }
  */
  /*
  template <typename T>
  bool has_overlap(
    S axa, S axb, S aya, S ayb, S aza, S azb,
    S bxa, S bxb, S bya, S byb, S bza, S bzb
  ) {
    if (axb <= bxa || bxb <= axa) return false;
    if (ayb <= bya || byb <= aya) return false;
    if (azb <= bza || bzb <= aza) return false;

    return true;
  }
  */
  template <typename T>
  bool has_overlap(
    S aax, S aay, S aaz, S abx, S aby, S abz,
    S bax, S bay, S baz, S bbx, S bby, S bbz
  ) {
    if (abx <= bax || bbx <= aax) return false;
    if (aby <= bay || bby <= aay) return false;
    if (abz <= baz || bbz <= aaz) return false;

    return true;
  }

  template <typename S>
  std::array<triple<S>, 2> free_cps_from_scratch(
    triple<std::vector<S>> &a, // The position of the box corner point.
    triple<std::vector<S>> &b, // The dimensions/sizes of the boxes.
    triple<S> c // The dimensions/sizes of the container.
  ) {
    // Create three fake boxes simulating the container walls.
    HBD_FOR3D(di) {
      HBD_FOR3D(dj) {
        if (di == dj) {
          a[di].push_back(0);
          b[di].push_back(1);
        } else {
          a[di].push_back(1);
          b[di].push_back(c[di]);
        }
      }
    }
    const auto free_cps = free_cps_from_scratch(a, b);
    HBD_FOR3D(d) {
      a[d].resize(a[d].size() - 3);
      b[d].resize(b[d].size() - 3);
    }

    return free_cps;
  }

  template <typename S>
  std::array<triple<std::vector<S>>, 2> free_cps_from_scratch(
    const triple<std::vector<S>> &a, // The position of the box corner point.
    const triple<std::vector<S>> &b  // The dimensions/sizes of the boxes.
  ) {
    // Change the names to make the code easier to read.
    const std::vector<S> &x = a[0], &y = a[1], &z = a[2];
    const std::vector<S> &p = b[0], &q = b[1], &r = b[2];
    const S l = c[0], w = c[1], h = c[2];

    assert( // All vectors have the same size.
      x.size() == y.size() && y.size() == z.size() && z.size() == p.size() &&
      p.size() == q.size() && q.size() == r.size()
    );

    const size_t n = x.size(); // number of boxes (i.e., size of all vectors)

    #ifndef NDEBUG
    // No overlap between boxes.
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = i + 1; j < n; ++j) {
        assert(!has_overlap(
          x[i], y[i], z[i], x[i] + p[i], y[i] + q[i], z[i] + r[i],
          x[j], y[j], z[j], x[j] + p[j], y[j] + q[j], z[j] + r[j]
        ));
      }
    }
    #endif
    
    std::array<triple<std::vector<S>>, 2> free_cps; // returned object
    for (size_t k = 0; k < n; ++k) {
      const S z_ = z[k] + h[k];
      for (size_t j = 0; j < n; ++j) {
        const S y_ = y[j] + w[j];
        if (z[j] + h[j] <= z_ || y[k] + w[k] <= y_) continue;
        for (size_t i = 0; i < n; ++i) {
          const S x_ = x[i] + l[i];
          if (
            x[k] + l[k] <= x_ || x[j] + l[j] <= x_ ||
            z[i] + h[i] <= z_ || y[i] + w[i] <= y_
          ) continue;
          const S p_ = max(x_, x[k], x[j]) - x_ + 1;
          const S q_ = max(y_, y[k], y[i]) - y_ + 1;
          const S r_ = max(z_, z[j], z[i]) - z_ + 1;
          for (size_t o = 0; o < n; ++o) {
            if (has_overlap(
              x_, y_, z_, x_ + p_, y_ + q_, z_ + r_,
              x[o], y[o], z[o], x[o] + l[o], y[o] + w[o], z[o] + h[o]
            )) goto after_pushing_cp;
          }
          free_cps[0][0].push_back(x_);
          free_cps[0][1].push_back(y_);
          free_cps[0][2].push_back(z_);
          free_cps[1][0].push_back(p_);
          free_cps[1][1].push_back(q_);
          free_cps[1][2].push_back(r_);
          after_pushing_cp:;
        }
      }
    }
    return free_cps;
  }

  // Common prefixes and suffixes:
  //  * c: container (preffix)
  //  * b: box (preffix)
  //  * cp: corner point (preffix)
  //  * nb: number (suffix)
  //
  // Template parameters:
  //  * an integer type large enough for the box/container sizes
  //  * an integer type large enough for the number of boxes available
  //  * an integer type large enough for the number of corner points created
  template <typename S, typename BI, typename CPI>
  struct space {
    // Container dimensions. Initialized with zeroes to make obvious that was
    // forgotten to be initialized.
    triple<S> c_dim{0, 0, 0};

    // Properties of the boxes.
    // The first position including the box.
    triple< std::vector<S> > b_start;
    // The first position after the box end.
    triple< std::vector<S> > b_after;
    // The indexes of the corner points clogged by the respective box.
    std::vector< std::vector<CPI> > b_clog_cp_ix;
    // The size of the corner-pointed-related vectors before the placement
    // of the respective box.
    std::vector< CPI > cp_nb_before_b;

    // Properties of the corner points.
    // The corner point coordinates.
    triple< std::vector<S> > cp;
    // Corner point minimal box after end. The first position after the
    // minimal box that can be fitted in the corner point.
    triple< std::vector<S> > cp_min_b_after;
    // Number of boxes that clog the corner point, making impossible to
    // place a box in the respective corner point. Note that a corner point
    // can be clogged just because there is a box already placed on it.
    std::vector<BI> clog_nb;

    space(triple<S> c_dim) :
      c_dim(c_dim),
      cp{{1}, {1}, {1}},
      cp_min_b_after{{1}, {1}, {1}},
      clog_nb{0}
      {}

    // auxiliar method
    void cp_resize(CPI size) {
      HBD_FOR3D(d) {
        cp[d].resize(size);
        cp_min_b_after[d].resize(size);
      }
      clog_nb.resize(size);
    }

    // auxiliar method
    void b_pop(void) {
      HBD_FOR3D(d) {
        b_start[d].pop();
        b_after[d].pop();
      }
      b_clog_cp_ix.pop();
      cp_nb_before_b.pop();
    }

    bool can_be_placed(triple<S> b_dim, CPI cp_ix) const {
      // Check if the cp is already known to be blocked by other boxes.
      if (clog_nb[cp_ix]) return false;

      HBD_FOR3D(d) {
        // Check if the box protrudes the container.
        if (cp[d][cp_ix] + b_dim[d] - 1 > c_dim[d]) return false;
        // Check if the box has the minimum size required by the cp.
        if (b_dim[d] < cp_min_b_after[d][cp_ix]) return false;
      }

      // Check if the box overlaps with other boxes (beyond the corner
      // point minimal box dimensions).
      const S aax = cp[0][cp_ix], aay = cp[1][cp_ix], aaz = cp[2][cp_ix],
        abx = aax + b_dim[0], aby = aay + b_dim[1], abz = aaz + b_dim[2];

      const BI b_nb = static_cast<BI>(b_start[0].size());
      for (BI i = 0; i < b_nb; ++i) {
        const S bax = b_start[0][i], bay = b_start[1][i], baz = b_start[2][i],
          bbx = b_after[0][i], bby = b_after[1][i], bbz = b_after[2][i];
        if (has_overlap(
          aax, aay, aaz, abx, aby, abz, bax, bay, baz, bbx, bby, bbz
        )) {
          return false;
        }
      }

      return true;
    }

    // Removes the last box from the structure (and any corner points
    // added by it, as also unclogging any corner points blocked by
    // the box existence). There is no subproduct of this computation
    // that is worth returning (everything relevant could be acessed
    // in O(1) before or after calling pop).
    void pop(void) {
      assert(cp_nb_before_b.size() > 0);
      cp_resize(cp_nb_before_b.back());
      const auto it_end = b_clog_cp_ix.back().end();
      auto it = b_clog_cp_ix.back().begin();
      while (it != it_end) --clog_nb[*it++];
      b_pop();
    }

    void place(triple<S> b_dim, CPI cp_ix) {
      assert(can_be_placed(b_dim, cp_ix));
      
      // Push the box into the space structure.
      HBD_FOR3D(d) {
        b_start[d].push_back(cp[d][cp_ix]);
        b_after[d].push_back(cp[d][cp_ix] + b_dim[d]);
      }
      // Saves the number of cps before the box was placed.
      cp_nb_before_b.push_back(static_cast<CPI>(cp[0].size()));
      // The new box starts with an empty list of clogged cps.
      b_clog_cp_ix.emplace_back();

      // Now we find which corner points the new box clogs.
      const S aax = cp[0][cp_ix], aay = cp[1][cp_ix], aaz = cp[2][cp_ix],
        abx = aax + b_dim[0], aby = aay + b_dim[1], abz = aaz + b_dim[2];
      const CPI cp_nb = static_cast<CPI>(cp[0].size());
      for (CPI i = 0; i < cp_nb; ++i) {
        const S bax = cp[0][i], bay = cp[1][i], baz = cp[0][i],
          bbx = cp_min_b_after[0][i], bby = cp_min_b_after[1][i],
          bbz = cp_min_b_after[2][i];
        if (has_overlap(
          aax, aay, aaz, abx, aby, abz, bax, bay, baz, bbx, bby, bbz
        )) {
          b_clog_cp_ix.back().push_back(i);
          ++clog_nb[i];
        }
      }
      
      // After placing a box in the cp, the cp should be marked as clogged.
      assert(clog_nb[cp_ix] > 0);

      // Create three fake boxes simulating the container walls.
      HBD_FOR3D(di) {
        HBD_FOR3D(dj) {
          if (di == dj) {
            b_start[di].push_back(0);
            b_after[di].push_back(1);
          } else {
            b_start[di].push_back(1);
            b_after[di].push_back(c_dim[di]);
          }
        }
      }
      
      // Number of boxes counting the newly placed and the three fake walls.
      const BI n = static_cast<BI>(b_start[0].size());

      // Change the names to make the code easier to read.
      const std::vector<S> &x = b_start[0], &y = b_start[1], &z = b_start[2];
      const std::vector<S> &x_ = b_after[0], &y_ = b_after[1],
        &z_ = b_after[2];
      HBD_FOR3D(d) {
        for (BI k = (d == 2 ? n - 4 : 0); k < (d == 2 ? n - 3 : n); ++k) {
          for (BI j = (d == 1 ? n - 4 : 0); j < (d == 1 ? n - 3 : n); ++j) {
            if (z_[j] <= z_[k] || y_[k] <= y_[j]) continue;
            for (BI i = (d == 0 ? n - 4 : 0); i < (d == 0 ? n - 3 : n); ++i) {
              if (
                x_[k] <= x_[i] || x_[j] <= x_[i] ||
                z_[i] <= z_[k] || y_[i] <= y_[j]
              ) continue;
              const S cx_ = max(x_[i], x[k], x[j]);
              const S cy_ = max(y_[j], y[i], y[k]);
              const S cz_ = max(z_[k], z[j], z[i]);
              for (BI o = 0; o < n - 4; ++o) {
                if (has_overlap(
                  x_[i], y_[j], z_[k], cx_, cy_, cz_,
                  x[o], y[o], z[o], x_[o], y_[o], z_[o]
                )) {
                  goto after_pushing_cp;
                }
              }
              cp[0].push_back(x_[i]);
              cp[1].push_back(y_[j]);
              cp[2].push_back(z_[k]);
              cp_min_b_after[0].push_back(cx_);
              cp_min_b_after[1].push_back(cy_);
              cp_min_b_after[2].push_back(cz_);
              // Every new cp starts unclogged.
              clog_nb.push_back(0);
              after_pushing_cp:;
            } // end of for x i
          } // end of for y j
        } // end of for z k
      } // end of HBD_FOR3D(d)
      // Remove the three fake boxes simulating container walls.
      HBD_FOR3D(d) {
        b_start[d].resize(n - 3);
        b_after[d].resize(n - 3);
      }
    } // end of method place
  }; // end of class space
} // end of namespace hbd

#endif // HBD_EXACT_CORNER_POINTS_HPP

