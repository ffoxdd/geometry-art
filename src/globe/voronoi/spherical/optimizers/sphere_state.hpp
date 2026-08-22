#ifndef GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_OPTIMIZERS_SPHERE_STATE_HPP_
#define GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_OPTIMIZERS_SPHERE_STATE_HPP_

#include "../../../types.hpp"
#include "../../../fields/spherical/region_integrals.hpp"
#include <cstddef>
#include <vector>

namespace globe::voronoi::spherical {

// One pass over the diagram: every cell's density integrals, and every
// bisector's, shared between the two cells that border it.
struct CellState {
    double mass;
    Vector3 first_moment;
};

struct EdgeState {
    size_t neighbor_index;
    fields::spherical::RegionIntegrals integrals;
};

struct SphereState {
    std::vector<CellState> cells;
    std::vector<std::vector<EdgeState>> edges;
};

} // namespace globe::voronoi::spherical

#endif //GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_OPTIMIZERS_SPHERE_STATE_HPP_
