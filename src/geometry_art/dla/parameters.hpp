#ifndef GEOMETRY_ART_DLA_PARAMETERS_HPP_
#define GEOMETRY_ART_DLA_PARAMETERS_HPP_

#include <cstddef>

namespace geometry_art::dla {

// overlap is a fraction of the particle radius; spawn_margin is measured in
// particle radii beyond the aggregate's reach.
struct Parameters {
    std::size_t particle_count = 2000;
    double particle_radius = 1.0;
    double overlap = 0.01;
    double spawn_margin = 2.0;

    [[nodiscard]] double overlap_distance() const {
        return overlap * particle_radius;
    }
};

} // namespace geometry_art::dla

#endif //GEOMETRY_ART_DLA_PARAMETERS_HPP_
