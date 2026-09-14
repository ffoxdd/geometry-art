#ifndef GEOMETRY_ART_DLA_CALLBACK_HPP_
#define GEOMETRY_ART_DLA_CALLBACK_HPP_

#include "aggregate.hpp"
#include "particle_index/particle_index.hpp"
#include "particle_index/rtree_particle_index.hpp"
#include <functional>

namespace geometry_art::dla {

template<ParticleIndex ParticleIndexType = RTreeParticleIndex>
using Callback = std::function<void(const Aggregate<ParticleIndexType>&)>;

template<ParticleIndex ParticleIndexType = RTreeParticleIndex>
[[nodiscard]] inline Callback<ParticleIndexType> noop_callback() {
    return [](const Aggregate<ParticleIndexType>&) {};
}

} // namespace geometry_art::dla

#endif //GEOMETRY_ART_DLA_CALLBACK_HPP_
