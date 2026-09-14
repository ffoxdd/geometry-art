#ifndef GEOMETRY_ART_TESTING_FLAT_SCATTER_HPP_
#define GEOMETRY_ART_TESTING_FLAT_SCATTER_HPP_

#include "../types.hpp"
#include "../geometry/planar/domain.hpp"
#include "../voronoi/flat/core/diagram.hpp"
#include <cmath>
#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

namespace geometry_art::testing {

using geometry::planar::Domain;

// A two-dimensional low-discrepancy scatter for flat tests. Rows of sites
// make nearly striped cells whose capacities are close to linear in the
// sites, which starves derivative tests of their subject; this spreads in
// both directions.
inline std::vector<Vector2> flat_scatter(size_t count, double width, double height) {
    std::vector<Vector2> sites;
    sites.reserve(count);

    for (size_t k = 0; k < count; ++k) {
        double x = std::fmod(0.13 + 0.7548776662466927 * static_cast<double>(k), 1.0) * width;
        double y = std::fmod(0.41 + 0.5698402909980532 * static_cast<double>(k), 1.0) * height;
        sites.emplace_back(x, y);
    }

    return sites;
}

inline std::unique_ptr<voronoi::flat::Diagram> diagram_of(std::vector<Vector2> sites, const Domain& domain) {
    return std::make_unique<voronoi::flat::Diagram>(domain, std::move(sites));
}

inline std::unique_ptr<voronoi::flat::Diagram> scattered_diagram(size_t count, const Domain& domain) {
    return diagram_of(flat_scatter(count, domain.width, domain.height), domain);
}

inline std::unique_ptr<voronoi::flat::Diagram> scattered_torus(size_t count, double width, double height) {
    return scattered_diagram(count, Domain::torus(width, height));
}

inline std::unique_ptr<voronoi::flat::Diagram> scattered_plane(size_t count, double width, double height) {
    return scattered_diagram(count, Domain::plane(width, height));
}

} // namespace geometry_art::testing

#endif //GEOMETRY_ART_TESTING_FLAT_SCATTER_HPP_
