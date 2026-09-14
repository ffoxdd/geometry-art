#ifndef GEOMETRY_ART_DLA_SPAWN_SHELL_KILLING_SPAWN_SHELL_HPP_
#define GEOMETRY_ART_DLA_SPAWN_SHELL_KILLING_SPAWN_SHELL_HPP_

#include "../direction_sampler/direction_sampler.hpp"
#include "../direction_sampler/uniform_direction_sampler.hpp"
#include "../../types.hpp"
#include <CGAL/assertions.h>
#include <utility>

namespace geometry_art::dla {

// The classic far-field rule, kept for comparison against the harmonic one:
// a walker wanders freely until it strays past kill_factor times the shell
// radius, then is killed and respawned. Approximate — a killed walker still
// had radius / distance odds of returning — with the bias shrinking as the
// kill factor grows and the walks lengthening with it.
template<DirectionSampler DirectionSamplerType = UniformDirectionSampler<>>
class KillingSpawnShell {
 public:
    explicit KillingSpawnShell(double kill_factor) :
        KillingSpawnShell(kill_factor, DirectionSamplerType()) {
    }

    KillingSpawnShell(double kill_factor, DirectionSamplerType direction_sampler) :
        _kill_factor(kill_factor),
        _direction_sampler(std::move(direction_sampler)) {
        CGAL_precondition(kill_factor >= 1.0);
    }

    [[nodiscard]] Vector3 spawn(double radius);
    [[nodiscard]] Vector3 confine(const Vector3& position, double radius);

 private:
    double _kill_factor;
    DirectionSamplerType _direction_sampler;
};

template<DirectionSampler DirectionSamplerType>
Vector3 KillingSpawnShell<DirectionSamplerType>::spawn(double radius) {
    return radius * _direction_sampler.sample();
}

template<DirectionSampler DirectionSamplerType>
Vector3 KillingSpawnShell<DirectionSamplerType>::confine(const Vector3& position, double radius) {
    if (position.norm() <= _kill_factor * radius) {
        return position;
    }

    return spawn(radius);
}

} // namespace geometry_art::dla

#endif //GEOMETRY_ART_DLA_SPAWN_SHELL_KILLING_SPAWN_SHELL_HPP_
