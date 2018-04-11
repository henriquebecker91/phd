#include <random>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <cassert>
#include <cstdbool>
#include <random>

#ifndef HBM_PRINT_VAR
  #define HBM_PRINT_VAR(var) out << #var ": " << var << std::endl
#endif

// For the classical DARP (multi-vehicle, homogeneous fleet, single depot)
// there can be more than one instance type or interesting data representation.
// Consequently, the solving method will not receive a concrete data structure
// but a templated type. Some things clearly can be stored as simple variables:
// number of requests (n), number of locations (nl), number of vehicles (m),
// vehicle capacity (vc), max route time (mrot). Such variables should always be
// present in the data structure. The time windows and service time will
// probably always be vectors, but they will have inline 'get' methods because
// the code can refer to locations using negative numbers or not (and sometimes
// all service times of an instance are equal, so they will not be a vector).
// Coordinates should not be directly referred in the code, as the problem
// respect the triangle unequality, but is not the metric DARP (some instances
// are metric, but there's no guarantee). The maximum ride time (mrit) is also
// behind a method because some definitions of the problem have different
// values for each request (so it would be a vector) and others use a single
// value.

// V for vehicle indexes
// L for location indexes
// P for people amounts and vehicle capacities
// C for costs between locations
// T for times between locations and time windows/amounts
// X for the value of XY coords
template <typename V, typename L, typename P, typename C, typename T, typename X>
struct CL2003_darp_inst_t {
  L n; /* number of requests */
  L nl; /* number of locations, always n*2 + 1 */
  V m; /* number of vehicles */
  P vc; /* vehicle capacity */
  T mrot; /* max route time */
  T mrit_; /* max ride time */

  std::vector<X> x;
  std::vector<X> y;
  std::vector<P> q_;
  std::vector<T> d_;
  std::vector<T> e_;
  std::vector<T> l_;

  /* ignores argument on purpose, all mrit is equal in these instances */
  inline T mrit(const L &a) { return mrit_; }

  inline T q(const L &a) { return q_[a]; }
  inline T d(const L &a) { return d_[a]; }
  inline T e(const L &a) { return e_[a]; }
  inline T l(const L &a) { return l_[a]; }

  inline C cost(const L &a, const L &b) {
    return static_cast<C>(
      std::sqrt(std::pow(x[a] - x[b], 2) + std::pow(y[a] - y[b], 2))
    );
  }

  inline T time(const L &a, const L &b) {
    return static_cast<T>(
      std::sqrt(std::pow(x[a] - x[b], 2) + std::pow(y[a] - y[b], 2))
    );
  }

  CL2003_darp_inst_t(std::istream &f) {
    f >> m >> n >> mrot >> vc >> mrit_;

    assert(n % 2 == 0);

    // n comes as the number of non-depot locations
    n = n/2;

    nl = 2*n + 1;

    x.reserve(nl);
    y.reserve(nl);
    q_.reserve(nl);
    d_.reserve(nl);
    e_.reserve(nl);
    l_.reserve(nl);

    L unused_id;
    X x_, y_;
    P q_in;
    T d_in, e_in, l_in;
    for (L i = 0; i < nl; ++i) {
      f >> unused_id >> x_ >> y_ >> q_in >> d_in >> e_in >> l_in;

      assert(unused_id == i);

      x.push_back(x_);
      y.push_back(y_);
      q_.push_back(q_in);
      d_.push_back(d_in);
      e_.push_back(e_in);
      l_.push_back(l_in);
    }
  }

  void debug_output(std::ostream &out = std::cout) {
    HBM_PRINT_VAR(n);
    HBM_PRINT_VAR(nl);
    HBM_PRINT_VAR(m);
    HBM_PRINT_VAR(vc);
    HBM_PRINT_VAR(mrot);
    HBM_PRINT_VAR(mrit_);
    for (L i = 0; i < nl; ++i) {
      out << "x[" << i << "]: " << x[i] << std::endl;
      out << "y[" << i << "]: " << y[i] << std::endl;
      out << "q(" << i << "): " << q(i) << std::endl;
      out << "d(" << i << "): " << d(i) << std::endl;
      out << "e(" << i << "): " << e(i) << std::endl;
      out << "l(" << i << "): " << l(i) << std::endl;
    }
  }
};

template <typename V, typename L, typename P, typename C, typename T>
using CL2003_darp_inst_wt = CL2003_darp_inst_t<V, L, P, C, T, float>;

template <typename K, typename V>
struct set_map_t;

template <
  template <typename, typename, typename, typename, typename> class O,
  typename V,
  typename L,
  typename P,
  typename C,
  typename T>
struct sol_t {
  std::vector<set_map_t<L, L>*> chosen_routes;
  O<V, L, P, C, T> const *inst;
  set_map_t<L, L>* all_routes;
  std::vector<V> route_with_req;
  std::vector<std::vector<L>> reqs_in_route;
  C cost;
  std::vector<C> req_cost;
  std::mt19937 mock_cost{0};
  std::uniform_int_distribution<int> mock_cost_dist{100, 10000};
  std::uniform_int_distribution<int> mock_feasibility_dist{0, 999999};
  V open_routes_var{1};

  sol_t(const O<V, L, P, C, T> * const inst_) : inst(inst_)
  {
    all_routes = new set_map_t<L, L>(inst->n);
    chosen_routes = std::vector<set_map_t<L, L> *>(inst->m, all_routes);
    route_with_req = std::vector<V>(inst->n + 1, inst->m);
    req_cost = std::vector<C>(inst->n + 1, 0);
    reqs_in_route = std::vector<std::vector<L>>(inst->m);
    for (auto &route : reqs_in_route) route.reserve(inst->n);
  }

  void remove_if_present(const L &r) {
    if (route_with_req[r] == inst->m) return;
    assert(reqs_in_route[route_with_req[r]].back() == r);
    cost -= req_cost[r];
    reqs_in_route[route_with_req[r]].pop_back();
    if (reqs_in_route[route_with_req[r]].empty()) --open_routes_var;
    route_with_req[r] = inst->m;
  }

  bool add_req_to_route(const L &req, const V &route) {
    assert(req <= inst->n);
    assert(route < inst->m);
    if (!mock_feasibility_dist(mock_cost)) return false;

    route_with_req[req] = route;
    if (reqs_in_route[route].empty() && open_routes_var < inst->m)
      ++open_routes_var;
    reqs_in_route[route].push_back(req);
    req_cost[req] = mock_cost_dist(mock_cost);
    cost += req_cost[req];

    return true;
  }

  void save_state(std::vector<set_map_t<L, L>*> &bkv_sol) {
    std::cout << "save_state:";
    for (auto &route : reqs_in_route) {
      for (auto &req : route) {
        std::cout << " " << std::to_string(req);
      }
      std::cout << " |";
    }
    std::cout << std::endl;
  }

  inline V open_routes(void) {
    return open_routes_var;
  }
};

template <
  template <typename, typename, typename, typename, typename> class O,
  typename V,
  typename L,
  typename P,
  typename C,
  typename T>
bool backtrack(
  L &curr_req,
  std::vector<V> &curr_try_route,
  sol_t<O, V, L, P, C, T> &curr_sol
) {
  // backtrack is only called after testing all distinct route insertions
  // of request number curr_route.size()
  #ifndef NDEBUG
  assert(curr_req > 0 && curr_req <= curr_sol.inst->n);
  if (curr_try_route[curr_req] != curr_sol.open_routes()) {
    auto &out = std::cout;
    HBM_PRINT_VAR(curr_try_route[curr_req]);
    HBM_PRINT_VAR(curr_sol.open_routes());
  }
  #endif
  curr_sol.remove_if_present(curr_req);
  curr_try_route[curr_req] = 0;
  while (true) {
    --curr_req;
    if (curr_req == 0) return false;
    if (curr_sol.reqs_in_route[curr_sol.route_with_req[curr_req]].size() == 1) continue;
    curr_sol.remove_if_present(curr_req);
    while (++curr_try_route[curr_req] < curr_sol.open_routes()) {
      if (curr_sol.add_req_to_route(curr_req, curr_try_route[curr_req]))
        return true;
    }
    if (curr_try_route[curr_req] >= curr_sol.open_routes()) {
      curr_try_route[curr_req] = 0;
    }
  }
}

template <
  template <typename, typename, typename, typename, typename> class O,
  typename V,
  typename L,
  typename P,
  typename C,
  typename T>
void test_all_mv_schedules(
  const O<V, L, P, C, T> &inst,
  std::vector<set_map_t<L, L>*> &bkv_sol,
  C &bkv
) {
  bool found_feasible_sol = false;
  L curr_req = 0;
  std::vector<V> curr_try_route(inst.n + 1, 0);
  sol_t<O, V, L, P, C, T> curr_sol(&inst);

  while (true) {
    if (curr_req < inst.n) {
      ++curr_req;
      for (; curr_try_route[curr_req] < curr_sol.open_routes(); ++curr_try_route[curr_req])
        if (curr_sol.add_req_to_route(curr_req, curr_try_route[curr_req]))
          break;
      if (curr_try_route[curr_req] == curr_sol.open_routes() &&
          !backtrack(curr_req, curr_try_route, curr_sol)) break;
    } else {
      assert(curr_req == inst.n);
      // the code only enters in this 'else' if the last loop increased
      // curr_route.size() to inst.m and it found a feasible solution before
      // curr_route.back() == open_routes
      if (!found_feasible_sol || bkv > curr_sol.cost) {
        bkv = curr_sol.cost;
        curr_sol.save_state(bkv_sol);
        found_feasible_sol = true;
      }
      curr_sol.remove_if_present(inst.n);
      const V open_routes = curr_sol.open_routes();
      for (; ++curr_try_route[curr_req] < open_routes; ) {
        if (curr_sol.add_req_to_route(curr_req, curr_try_route[curr_req]) &&
            (!found_feasible_sol || bkv > curr_sol.cost)) {
          bkv = curr_sol.cost;
          found_feasible_sol = true;
        }
        curr_sol.save_state(bkv_sol);
        curr_sol.remove_if_present(curr_req);
      }
      if (!backtrack(curr_req, curr_try_route, curr_sol)) break;
    }
  }
}

template <typename K, typename V>
struct set_map_t {
  std::vector<set_map_t<K, V>*> c;
  V v;

  set_map_t(
    const size_t &set_size,
    set_map_t<K, V>* const &parent = nullptr,
    const V &v_ = V()) :
    c(set_size, nullptr), v(v_)
  {
    assert(set_size > 0);
    c[0] = parent;
  }

  template <typename I>
  set_map_t<K, V>* walk_existing_path(I &begin, const I &end) const {
    set_map_t<K, V>* curr = this;
    while (curr.c[*begin] != nullptr && begin != end) {
      curr = curr.c[*(begin++)];
    }
    return curr;
  }

  template <typename I>
  bool exists(I &begin, const I &end) {
    // walk_existing_path is used below only to advance begin
    static_cast<void>(walk_existing_path(begin, end));
    return begin == end;
  }

  template <typename I>
  set_map_t<K, V>* forge_path(I &begin, const I &end) {
    set_map_t<K, V>* curr = walk_existing_path(begin, end);
    // creates the remainder of the path if necessary
    for (; begin != end; ++begin) {
      curr.c[*begin] = new set_map_t<K, V>(*begin, curr);
      curr = curr.c[*begin];
    }
    return curr;
  }

  template <typename I>
  V& get(I &begin, const I &end) {
    auto node = forge_path(begin, end);

    return node->v;
  }

  V& operator[](const std::vector<K> &set) {
    return get(set.begin(), set.end());
  }
};

// OLD CODE THAT MAY OR MAY NOT BE USED
/*
template <typename X>
struct coord2D_t {
  X x{0};
  X y{0};

  X dist_to(const coord2D_t<X> &o) const {
    return static_cast<X>(sqrt(pow(x - o.x, 2) + pow(y - o.y, 2)));
  }

  coord2D_t(void) : x(0), y(0) {}
  coord2D_t(const X &x_, const X &y_) : x(x_), y(y_) {}
};


template <typename W, typename P>
struct request_t {
  coord2D_t<W> pickup;
  coord2D_t<W> delivery;
  P qt_people;

  request_t(
    const coord2D_t<W> pickup_,
    const coord2D_t<W> delivery_,
    const P qt_people_) :
    pickup(pickup_),
    delivery(delivery_),
    qt_people(qt_people_) {}
};

template <typename W, typename P>
struct vehicle_t {
  coord2D_t<W> start;
  coord2D_t<W> end;
  P capacity;

  vehicle_t(
    const coord2D_t<W> start_,
    const coord2D_t<W> end_,
    const P capacity_) :
    start(start_),
    end(end_),
    capacity(capacity_) {}
};

// Things needed to create a simple random coord instance:
// max_x, max_y, num_vehicles, vehicle_coords, num_requests,
// request_{pickup_point,delivery_point,qt_people}
template <typename V, typename L, typename W, typename P>
struct coord_darp_inst_t : darp_inst_t<V, L, W> {
  std::vector< vehicle_t<W, P> > vehicles;
  std::vector< request_t<W, P> > requests;
  // https://isocpp.org/wiki/faq/virtual-functions#virtual-ctors
  virtual darp_inst_t<V, L, W>* clone(void) const {
    return nullptr; // not implemented yet
  };

  virtual ~darp_inst_t(std::istream in) {
    
  };

  // RNG -- a random number engine
  // DVC -- a distribution of vehicle capacities
  // DRP -- a distribution of people in requests
  template<typename RNG, typename DVC, typename DPR>
  coord_darp_inst_t(
    const V num_vehicles,
    const L num_requests,
    const W max_x,
    const W max_y,
    RNG &rng,
    DVC &dvc,
    DPR &dpr
  ) {
    std::uniform_int_distribution<W> dx(0, max_x);
    std::uniform_int_distribution<W> dy(0, max_y);
    for (V i = 0; i < num_vehicles; ++i) {
      vehicles.emplace_back(
        coord2D_t<W>(dx(rng), dy(rng)),
        coord2D_t<W>(dx(rng), dy(rng)),
        dvc(rng)
      );
    }
    for (L i = 0; i < num_requests; ++i) {
      requests.emplace_back(
        coord2D_t<W>(dx(rng), dy(rng)),
        coord2D_t<W>(dx(rng), dy(rng)),
        dpr(rng)
      );
    }
  }

  std::string to_dot(void) {
    std::string s = "digraph G {\n";
    
    for (V i = 1; i <= vehicles.size(); ++i) {
      s += "node_" + std::to_string(i) + 
           "[style=filled,fillcolor=chartreuse, label=\"" +
           std::to_string(i) + "\", shape=circle,pos=\"" +
           std::to_string(vehicles[i-1].start.x) + "," +
           std::to_string(vehicles[i-1].start.y) + "!\"]\n";
      s += "node_" + std::to_string(vehicles.size() + i) + 
           "[style=filled,fillcolor=firebrick1, label=\"" +
           std::to_string(i) + "\", shape=circle,pos=\"" +
           std::to_string(vehicles[i-1].end.x) + "," +
           std::to_string(vehicles[i-1].end.y) + "!\"]\n";
      s += "node_" + std::to_string(i) + " -> node_" +
           std::to_string(vehicles.size() + i) + "\n";
    }
    const size_t offset = 2*vehicles.size();
    for (L i = 1; i <= requests.size(); ++i) {
      s += "node_" + std::to_string(offset + i) + 
           "[style=filled,fillcolor=white, label=\"" +
           std::to_string(i) + "\", shape=circle,pos=\"" +
           std::to_string(requests[i-1].pickup.x) + "," +
           std::to_string(requests[i-1].pickup.y) + "!\"]\n";
      s += "node_" + std::to_string(offset + requests.size() + i) + 
           "[style=filled,fillcolor=gray, label=\"" +
           std::to_string(i) + "\", shape=circle,pos=\"" +
           std::to_string(requests[i-1].delivery.x) + "," +
           std::to_string(requests[i-1].delivery.y) + "!\"]\n";
      s += "node_" + std::to_string(offset + i) + " -> node_" +
           std::to_string(offset + requests.size() + i) + "\n";
    }

    s += "}\n";

    return s;
  }

  std::string to_s(void) {
    std::string s = std::to_string(requests.size()) + "\t" +
                    std::to_string(vehicles.size()) + "\n";
    for (auto &v : vehicles) {
      s += std::to_string(v.capacity) + "\t" +
           std::to_string(v.start.x) + "\t" +
           std::to_string(v.start.y) + "\t" +
           std::to_string(v.end.x) + "\t" +
           std::to_string(v.end.y) + "\n";
    }
    for (auto &r : requests) {
      s += std::to_string(r.qt_people) + "\t" +
           std::to_string(r.pickup.x) + "\t" +
           std::to_string(r.pickup.y) + "\t" +
           std::to_string(r.delivery.x) + "\t" +
           std::to_string(r.delivery.y) + "\n";
    }
    return s;
  }
};

// V for vehicle indexes
// L for location indexes
// W for weights between locations
// P for people amounts
template <typename V, typename L, typename W, typename P>
struct darp_inst_t {
  //virtual ~darp_inst_t() {};
  //virtual darp_inst_t& operator=(darp_inst_t &b) = 0;

  // https://isocpp.org/wiki/faq/virtual-functions#virtual-ctors
  virtual darp_inst_t<V, L, W>* clone(void) const = 0;
  //virtual darp_inst_t* create() const = 0; what arguments it would take?

  //virtual W distance(const L a, const L b) const;
};
*/

