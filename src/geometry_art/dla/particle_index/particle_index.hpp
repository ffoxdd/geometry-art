#ifndef GEOMETRY_ART_DLA_PARTICLE_INDEX_PARTICLE_INDEX_HPP_
#define GEOMETRY_ART_DLA_PARTICLE_INDEX_PARTICLE_INDEX_HPP_

#include "../../types.hpp"
#include <concepts>
#include <cstddef>

namespace geometry_art::dla {

struct NearestParticle {
    std::size_t index;
    double distance;
};

// A dynamic nearest-neighbor index over particle centers: insertion one
// particle at a time as the aggregate grows, nearest queries between
// insertions.
template<typename T>
concept ParticleIndex = requires(T index, const Vector3& point, std::size_t particle) {
    { index.insert(point, particle) };
    { index.nearest(point) } -> std::convertible_to<NearestParticle>;
};

} // namespace geometry_art::dla

#endif //GEOMETRY_ART_DLA_PARTICLE_INDEX_PARTICLE_INDEX_HPP_
