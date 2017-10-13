#include <random>
#include <iostream>
#include <vector>
#include <string>

// V for vehicle indexes
// L for location indexes
// W for weights between locations
// P for people amounts
template <typename V, typename L, typename W>
struct darp_inst_t {
  //virtual ~darp_inst_t() {};
  //virtual darp_inst_t& operator=(darp_inst_t &b) = 0;

  // https://isocpp.org/wiki/faq/virtual-functions#virtual-ctors
  virtual darp_inst_t<V, L, W>* clone(void) const = 0;
  //virtual darp_inst_t* create() const = 0; what arguments it would take?

  //virtual W distance(const L a, const L b) const;
};

template <typename W>
struct coord2D_t {
  W x{0};
  W y{0};

  W distance_to(const coord2D_t<W> &o) const {
    return static_cast<W>(sqrt(pow(x - o.x, 2) + pow(y - o.y, 2)));
  }

  coord2D_t(const W x_, const W y_) : x(x_), y(y_) {}
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

  /*virtual ~darp_inst_t(std::istream in) {
    
  };*/

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

