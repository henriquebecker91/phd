#ifndef HBM_CHEN1995_MODEL_HPP
#define HBM_CHEN1995_MODEL_HPP

#include <cassert>
#include <cmath>
#include <ilcplex/ilocplex.h>
#include "ukp_common.hpp"

#ifndef HBM_PRINT_VAR
  #define HBM_PRINT_VAR(var) out << (#var ": ") << var << std::endl
#endif

namespace hbm {
  // Q: quantities, as the number of boxes (n) or containers (m)
  // C: coordinates, positions in an axis or distance in 1D~3D space
  template<typename Q, typename C>
  void chen1995(
    // Number of boxes (size of s, p, q, r, lx/y/z, wx/y/z, hx/y/z, a to f).
    const Q N,
    // Number of containers (size of the arrays inside s, size of n).
    const Q m,
    // Big-m (M) will be computed inside the procedure.
    // const C M,
    // Box-Container flag. One if box i is in container j; otherwise zero.
    bool **s,
    // Container-Use flag. One if container j has a box inside; otherise zero.
    bool *n,
    // lenght (p), width (q), height (r) of the boxes
    const C const *p, const C const *q, const C const *r,
    // lenght (L), width (W), height (H) of the containers
    const C const *L, const C const *W, const C const *H,
    // The position assigned to each box by the model.
    C *x, C *y, C *z,
    // The rotation assigned to each box by the model.
    // ex.: if wy[i] == 1 then the width of item i is aligned with the y axis
    C *lx, C *ly, C *lz, C *wx, C *wy, C *wz, C *hx, C *hy, C *hz,
    // The position of a box relative to other, as assigned by the model.
    // If a[i][k] == 1, then box i is at the left of box k, the letters meaning
    // is: left (a), right (b), behind (c), front (d), below (e), above (f).
    bool **a, bool **b, bool **c, bool **d, bool **d, bool **e, bool **f
  ) {
    using namespace std;
    //auto &out = cout;

    IloEnv env;

    //auto toIloIntArray [&env,&N] () { }

    //IloInt abv = 0; // all boxes volume
    //for (IloInt i = 0; i < N; ++i) abv += static_cast<IloInt>(p[i]*q[i]*r[i]);

    // container volumes
    IloIntArray cv(env, m);
    for (IloInt i = 0; i < m; ++i) cv[i] = static_cast<IloInt>(L[i]*W[i]*H[i]);

    IloNumVarArray n_(env, m, 0, 1, ILOINT);
    IloArray<IloNumVarArray> s_(env, N);
    for (IloInt i = 0; i < m; ++i) cv[i] = static_cast<IloInt>(L[i]*W[i]*H[i]);

    IloNumVarArray s_(env, m, 0, 1, ILOINT);
    for (IloInt i = 0; i < N; ++i) s_[i] = IloIntVarArray(env, m, 0, 1);

    // create MAX_X, MAX_Y and MAX_Z with the bound less largest box lenght

    IloNumVarArray x_(env, N, 0, MAX_X);
    IloNumVarArray y_(env, N, 0, MAX_Y);
    IloNumVarArray z_(env, N, 0, MAX_Z);

    IloIntVarArray lx_(env, N, 0, 1);
    IloIntVarArray ly_(env, N, 0, 1);
    IloIntVarArray lz_(env, N, 0, 1);
    IloIntVarArray wx_(env, N, 0, 1);
    IloIntVarArray wy_(env, N, 0, 1);
    IloIntVarArray wz_(env, N, 0, 1);
    IloIntVarArray hx_(env, N, 0, 1);
    IloIntVarArray hy_(env, N, 0, 1);
    IloIntVarArray hz_(env, N, 0, 1);
    // Only for reference of what needs to become a variable.
    //C *lx, C *ly, C *lz, C *wx, C *wy, C *wz, C *hx, C *hy, C *hz,
    //bool **a, bool **b, bool **c, bool **d, bool **d, bool **e, bool **f

    IloModel model(env);
    //cout << "before maximize" << endl;

    // In the paper, the objective function has a constant term (the sum of the
    // volume of all boxes) which is not really necessary for the model.
    model.add(IloMinimize(env, IloScalProd(cv, n_)));
    //model.add(IloMinimize(env, IloScalProd(cv, n_) - abv));

    for (IloInt k = 1; k < N; ++k) {
      for (IloInt i = 0; i < k; ++i) {
        model.add(xi + );
      }
    }

    IloCplex cplex(model);

    // configure cplex solver
    //cplex.setOut(env.getNullStream()); // disable output
    cplex.setParam(IloCplex::Param::RandomSeed, 0);
    // The AbsMIPGap is the absolute gap between best solution found and
    // the optimistic guess (upper or lower bound). The default is 1e-06
    // and we prefer to not change it because it can mess the CPLEX performance
    // without giving any precision gain. The MIPGap is a relative value
    // with default 1e-04, so any model with an objective value over 0.01
    // would stop first because of this relative gape than because of
    // the absolute one.
    //cplex.setParam(IloCplex::Param::MIP::Tolerances::AbsMIPGap, 0.0);
    cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, 0.0);

    // The cplex guarantees that setting the threads to one gives a
    // completely deterministic execution, so the Parallel param does
    // not need to be used to guarantee determinism.
    cplex.setParam(IloCplex::Param::Threads, 1);
    //cplex.setParam(IloCplex::Param::Parallel, 1);

    // Maybe define an internal time limit?
    cplex.setParam(IloCplex::Param::TimeLimit, 1800);
    //cplex.setParam(IloCplex::Param::MIP::Limits::TreeMemory, 7*1024);
    // Guarantee CPLEX is using the wall clock time (default).
    cplex.setParam(IloCplex::Param::ClockType, 2);
    // TODO: make some preliminary tests with the parameters below to decide
    // if one of them should be used.
    cplex.setParam(IloCplex::Param::Emphasis::MIP, 1);
    // Values for IloCplex::Param::Emphasis::MIP
    // 0 balance optimality and feasibility
    // 1 feasibility over optimality
    // 2 optimality over feasibility
    // 3 "even greater emphasis is placed on proving optimality"
    // 4 "consider this setting when the FEASIBILITY setting has
    //    difficulty finding solutions of acceptable quality."

    // Granularity of the information display by CPLEX.
    cplex.setParam(IloCplex::Param::Simplex::Display, 2);
    cplex.setParam(IloCplex::Param::MIP::Tolerances::Integrality, 0.0);

    // TODO: check if "numerical precision emphasis" should be set to one,
    // the default is zero that is to not worry much about numerical precision.
    // TODO: "presolve switch" maybe there is no reason to try to simplify
    // the model
    // TODO: check if "CPU mask to bind threads to cores" should be used
    // instead of taskset
    // TODO: check if CPLEX "integrality tolerance" has to be configured too.

    //cout << "before solve" << endl;
    cplex.solve();

    // The variables are from an IloIntVarArray, and even so their values need
    // to be captured in IloNumArray variable (i.e., a floating point variable)
    // and IloRound be used to get the real integer values.
    IloNumArray xv(env, n);
    cplex.getValues(xv, x);

    sol.opt = static_cast<P>(IloRound(cplex.getObjValue()));
    for (int i = 0; i < n; ++i) {
      if (IloRound(xv[i]) >= 1) sol.used_items.emplace_back(
        ukpi.items[i], static_cast<W>(IloRound(xv[i])), i
      );
    }

    env.end();
  }
}

#endif //HBM_CPLEX_UKP_MODEL_HPP

