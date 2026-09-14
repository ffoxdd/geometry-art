#ifndef GEOMETRY_ART_IO_SNAPSHOT_AGGREGATE_SNAPSHOT_HPP_
#define GEOMETRY_ART_IO_SNAPSHOT_AGGREGATE_SNAPSHOT_HPP_

#include "../../dla/aggregate.hpp"
#include "../../dla/particle.hpp"
#include "../../dla/particle_index/particle_index.hpp"
#include <vector>

namespace geometry_art::io::snapshot {

// An aggregate as plain data, for the writers and for polling viewers.
struct AggregateSnapshot {
    double particle_radius = 0.0;
    double reach = 0.0;
    std::vector<dla::Particle> particles;
};

template<dla::ParticleIndex ParticleIndexType>
[[nodiscard]] AggregateSnapshot capture(const dla::Aggregate<ParticleIndexType>& aggregate) {
    return AggregateSnapshot{
        aggregate.particle_radius(),
        aggregate.reach(),
        aggregate.particles()
    };
}

} // namespace geometry_art::io::snapshot

#endif //GEOMETRY_ART_IO_SNAPSHOT_AGGREGATE_SNAPSHOT_HPP_
