#ifndef GLOBEART_SRC_GLOBE_TESTING_FLAT_SCATTER_HPP_
#define GLOBEART_SRC_GLOBE_TESTING_FLAT_SCATTER_HPP_

#include "../types.hpp"
#include "../voronoi/flat/core/torus.hpp"
#include <cmath>
#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

namespace globe::testing {

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

inline std::unique_ptr<voronoi::flat::Torus> torus_of(
    std::vector<Vector2> sites,
    double width,
    double height
) {
    return std::make_unique<voronoi::flat::Torus>(width, height, std::move(sites));
}

inline std::unique_ptr<voronoi::flat::Torus> scattered_torus(size_t count, double width, double height) {
    return torus_of(flat_scatter(count, width, height), width, height);
}

} // namespace globe::testing

#endif //GLOBEART_SRC_GLOBE_TESTING_FLAT_SCATTER_HPP_
