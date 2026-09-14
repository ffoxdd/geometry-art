#ifndef GEOMETRY_ART_DLA_PARTICLE_HPP_
#define GEOMETRY_ART_DLA_PARTICLE_HPP_

#include "../types.hpp"
#include <cstddef>
#include <optional>

namespace geometry_art::dla {

// A frozen sphere of the aggregate: its center, and the index of the particle
// it stuck to. Only the seed particle has no parent.
struct Particle {
    Vector3 center;
    std::optional<std::size_t> parent;
};

} // namespace geometry_art::dla

#endif //GEOMETRY_ART_DLA_PARTICLE_HPP_
