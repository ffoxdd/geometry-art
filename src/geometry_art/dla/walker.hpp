#ifndef GEOMETRY_ART_DLA_WALKER_HPP_
#define GEOMETRY_ART_DLA_WALKER_HPP_

#include "aggregate.hpp"
#include "parameters.hpp"
#include "particle.hpp"
#include "particle_index/particle_index.hpp"
#include "direction_sampler/direction_sampler.hpp"
#include "direction_sampler/uniform_direction_sampler.hpp"
#include "spawn_shell/spawn_shell.hpp"
#include "spawn_shell/harmonic_spawn_shell.hpp"
#include "../types.hpp"
#include <CGAL/assertions.h>
#include <utility>

namespace geometry_art::dla {

// One particle's journey from the spawn shell to its place in the aggregate:
// jump by the clearance (exact walk-on-spheres), resolve any excursion beyond
// the shell, stick on reaching the overlap shell around the aggregate.
template<
    DirectionSampler DirectionSamplerType = UniformDirectionSampler<>,
    SpawnShell SpawnShellType = HarmonicSpawnShell<>
>
class Walker {
 public:
    explicit Walker(const Parameters& parameters) :
        Walker(parameters, DirectionSamplerType(), SpawnShellType()) {
    }

    Walker(const Parameters& parameters, DirectionSamplerType direction_sampler, SpawnShellType spawn_shell) :
        _parameters(parameters),
        _direction_sampler(std::move(direction_sampler)),
        _spawn_shell(std::move(spawn_shell)) {
        CGAL_precondition(parameters.overlap_distance() > 0.0);
    }

    template<ParticleIndex ParticleIndexType>
    [[nodiscard]] Particle settle(const Aggregate<ParticleIndexType>& aggregate);

 private:
    Parameters _parameters;
    DirectionSamplerType _direction_sampler;
    SpawnShellType _spawn_shell;

    template<ParticleIndex ParticleIndexType>
    [[nodiscard]] double spawn_radius(const Aggregate<ParticleIndexType>& aggregate) const;

    [[nodiscard]] Vector3 jump(const Vector3& position, double clearance);

    template<ParticleIndex ParticleIndexType>
    [[nodiscard]] Particle snapped(const Aggregate<ParticleIndexType>& aggregate, const Vector3& position) const;
};

template<DirectionSampler DirectionSamplerType, SpawnShell SpawnShellType>
template<ParticleIndex ParticleIndexType>
Particle Walker<DirectionSamplerType, SpawnShellType>::settle(const Aggregate<ParticleIndexType>& aggregate) {
    double shell_radius = spawn_radius(aggregate);
    Vector3 position = _spawn_shell.spawn(shell_radius);

    while (true) {
        double clearance = aggregate.clearance(position);

        if (clearance <= _parameters.overlap_distance()) {
            return snapped(aggregate, position);
        }

        position = _spawn_shell.confine(jump(position, clearance), shell_radius);
    }
}

template<DirectionSampler DirectionSamplerType, SpawnShell SpawnShellType>
template<ParticleIndex ParticleIndexType>
double Walker<DirectionSamplerType, SpawnShellType>::spawn_radius(
    const Aggregate<ParticleIndexType>& aggregate
) const {
    return aggregate.reach() + _parameters.spawn_margin * _parameters.particle_radius;
}

// A jump of exactly the clearance is exact, not heuristic: Brownian motion
// exits the largest empty ball around its start uniformly.
template<DirectionSampler DirectionSamplerType, SpawnShell SpawnShellType>
Vector3 Walker<DirectionSamplerType, SpawnShellType>::jump(const Vector3& position, double clearance) {
    return position + clearance * _direction_sampler.sample();
}

// The final placement: onto the contact sphere of the particle it hit, pushed
// in by the overlap distance, along the ray from that particle through the
// walker.
template<DirectionSampler DirectionSamplerType, SpawnShell SpawnShellType>
template<ParticleIndex ParticleIndexType>
Particle Walker<DirectionSamplerType, SpawnShellType>::snapped(
    const Aggregate<ParticleIndexType>& aggregate,
    const Vector3& position
) const {
    NearestParticle nearest = aggregate.nearest_particle(position);
    const Vector3& contact_center = aggregate.particles()[nearest.index].center;

    Vector3 direction = (position - contact_center).normalized();
    double stuck_distance = aggregate.contact_distance() - _parameters.overlap_distance();

    return Particle{contact_center + stuck_distance * direction, nearest.index};
}

} // namespace geometry_art::dla

#endif //GEOMETRY_ART_DLA_WALKER_HPP_
