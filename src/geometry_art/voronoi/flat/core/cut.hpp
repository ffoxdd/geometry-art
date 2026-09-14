#ifndef GEOMETRY_ART_VORONOI_FLAT_CORE_CUT_HPP_
#define GEOMETRY_ART_VORONOI_FLAT_CORE_CUT_HPP_

#include "../../../types.hpp"
#include <cstddef>
#include <variant>

namespace geometry_art::voronoi::flat {

// The two things that can bound a cell on a flat domain: the bisector with
// a neighbouring site, placed in the cell's own chart, or a wall of the
// domain. A cell edge lies on one cut; a cell vertex is where two meet, so
// the cut on the far side of a vertex is what pins it besides the edge's
// own bisector.
struct Bisector {
    size_t neighbor_index;
    Vector2 neighbor_position;
};

struct Wall {
    Vector2 inward_normal;
};

using Cut = std::variant<Bisector, Wall>;

} // namespace geometry_art::voronoi::flat

#endif //GEOMETRY_ART_VORONOI_FLAT_CORE_CUT_HPP_
