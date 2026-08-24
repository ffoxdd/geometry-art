#ifndef GLOBEART_SRC_GLOBE_VORONOI_STATE_HPP_
#define GLOBEART_SRC_GLOBE_VORONOI_STATE_HPP_

#include "../types.hpp"
#include "../fields/region_integrals.hpp"
#include <cstddef>
#include <vector>

namespace globe::voronoi {

// One pass over the diagram: every cell's density integrals, and every
// bisector's, shared between the two cells that border it. A cell's moments
// are taken in its own chart -- on a wrapped domain a cell is placed around
// its site -- so they are only ever paired with that site.
struct CellState {
    double mass;
    Vector3 first_moment;

    // The integral of rho |x|^2, which the CVT energy needs; on the unit
    // sphere it equals the mass.
    double squared_norm_moment;
};

// Every stored quantity is invariant under a common translation of the
// bisector and its two sites, so cells on a wrapped domain read the same
// values from charts a period apart.
struct EdgeState {
    size_t neighbor_index;
    double separation;
    double mass;
    Vector3 moment_about_own;
    Vector3 moment_about_neighbor;
};

struct DiagramState {
    std::vector<CellState> cells;
    std::vector<std::vector<EdgeState>> edges;
};

} // namespace globe::voronoi

#endif //GLOBEART_SRC_GLOBE_VORONOI_STATE_HPP_
