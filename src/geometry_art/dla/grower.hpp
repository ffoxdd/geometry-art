#ifndef GEOMETRY_ART_DLA_GROWER_HPP_
#define GEOMETRY_ART_DLA_GROWER_HPP_

#include "aggregate.hpp"
#include "callback.hpp"
#include "parameters.hpp"
#include "walker.hpp"
#include "particle_index/particle_index.hpp"
#include "particle_index/rtree_particle_index.hpp"
#include "direction_sampler/direction_sampler.hpp"
#include "direction_sampler/uniform_direction_sampler.hpp"
#include "spawn_shell/spawn_shell.hpp"
#include "spawn_shell/harmonic_spawn_shell.hpp"
#include <utility>

namespace geometry_art::dla {

// The outer loop: settle one walker at a time onto the aggregate until it
// holds the requested number of particles, reporting after each freeze.
template<
    ParticleIndex ParticleIndexType = RTreeParticleIndex,
    DirectionSampler DirectionSamplerType = UniformDirectionSampler<>,
    SpawnShell SpawnShellType = HarmonicSpawnShell<>
>
class Grower {
 public:
    Grower(
        const Parameters& parameters,
        Aggregate<ParticleIndexType> aggregate,
        Walker<DirectionSamplerType, SpawnShellType> walker,
        Callback<ParticleIndexType> callback
    ) :
        _parameters(parameters),
        _aggregate(std::move(aggregate)),
        _walker(std::move(walker)),
        _callback(std::move(callback)) {
    }

    [[nodiscard]] Aggregate<ParticleIndexType> grow();

 private:
    Parameters _parameters;
    Aggregate<ParticleIndexType> _aggregate;
    Walker<DirectionSamplerType, SpawnShellType> _walker;
    Callback<ParticleIndexType> _callback;
};

template<ParticleIndex ParticleIndexType, DirectionSampler DirectionSamplerType, SpawnShell SpawnShellType>
Aggregate<ParticleIndexType> Grower<ParticleIndexType, DirectionSamplerType, SpawnShellType>::grow() {
    while (_aggregate.size() < _parameters.particle_count) {
        _aggregate.freeze(_walker.settle(_aggregate));
        _callback(_aggregate);
    }

    return std::move(_aggregate);
}

} // namespace geometry_art::dla

#endif //GEOMETRY_ART_DLA_GROWER_HPP_
