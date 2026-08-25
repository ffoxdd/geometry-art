#ifndef GEOMETRY_ART_VORONOI_SPHERICAL_CORE_DIAGRAM_HPP_
#define GEOMETRY_ART_VORONOI_SPHERICAL_CORE_DIAGRAM_HPP_

#include "sphere.hpp"
#include <concepts>
#include <cstddef>
#include <vector>

namespace geometry_art::voronoi::spherical {

// What the integration and derivative machinery reads from a tessellation.
// Everything above this line is arithmetic on cells and bisectors, so a
// hand-built stand-in with a few known cells can drive it without a
// triangulation.
template<typename T>
concept Diagram = requires(const T& diagram, size_t index) {
    { diagram.size() } -> std::convertible_to<size_t>;
    { diagram.site(index) } -> std::convertible_to<cgal::Point3>;
    { diagram.cell(index) } -> std::convertible_to<Polygon>;
    { diagram.cell_edges(index) } -> std::convertible_to<std::vector<CellEdgeInfo>>;
};

static_assert(Diagram<Sphere>);

} // namespace geometry_art::voronoi::spherical

#endif //GEOMETRY_ART_VORONOI_SPHERICAL_CORE_DIAGRAM_HPP_
