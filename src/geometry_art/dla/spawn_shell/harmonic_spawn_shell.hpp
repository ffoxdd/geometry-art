#ifndef GEOMETRY_ART_DLA_SPAWN_SHELL_HARMONIC_SPAWN_SHELL_HPP_
#define GEOMETRY_ART_DLA_SPAWN_SHELL_HARMONIC_SPAWN_SHELL_HPP_

#include "harmonic_hit.hpp"
#include "../direction_sampler/direction_sampler.hpp"
#include "../direction_sampler/uniform_direction_sampler.hpp"
#include "../../types.hpp"
#include "../../math/interval.hpp"
#include "../../math/interval_sampler/interval_sampler.hpp"
#include "../../math/interval_sampler/uniform_interval_sampler.hpp"
#include <Eigen/Geometry>
#include <cmath>
#include <utility>

namespace geometry_art::dla {

// The exact far-field rule: a walker beyond the shell has its whole excursion
// resolved in one draw. It escapes to infinity with the classical probability
// and respawns fresh, or it returns and lands where the exterior Poisson
// kernel says a returning Brownian walker lands. No kill radius, no bias.
template<
    DirectionSampler DirectionSamplerType = UniformDirectionSampler<>,
    IntervalSampler IntervalSamplerType = UniformIntervalSampler
>
class HarmonicSpawnShell {
 public:
    HarmonicSpawnShell() = default;

    HarmonicSpawnShell(DirectionSamplerType direction_sampler, IntervalSamplerType interval_sampler) :
        _direction_sampler(std::move(direction_sampler)),
        _interval_sampler(std::move(interval_sampler)) {
    }

    [[nodiscard]] Vector3 spawn(double radius);
    [[nodiscard]] Vector3 confine(const Vector3& position, double radius);

 private:
    DirectionSamplerType _direction_sampler;
    IntervalSamplerType _interval_sampler;

    [[nodiscard]] bool returns(double distance, double radius);
    [[nodiscard]] Vector3 landing(const Vector3& position, double distance, double radius);
};

// Walkers from infinity hit an enclosing sphere uniformly, so entering and
// escaping-then-reentering draw from the same distribution.
template<DirectionSampler DirectionSamplerType, IntervalSampler IntervalSamplerType>
Vector3 HarmonicSpawnShell<DirectionSamplerType, IntervalSamplerType>::spawn(double radius) {
    return radius * _direction_sampler.sample();
}

template<DirectionSampler DirectionSamplerType, IntervalSampler IntervalSamplerType>
Vector3 HarmonicSpawnShell<DirectionSamplerType, IntervalSamplerType>::confine(
    const Vector3& position,
    double radius
) {
    double distance = position.norm();

    if (distance <= radius) {
        return position;
    }

    if (!returns(distance, radius)) {
        return spawn(radius);
    }

    return landing(position, distance, radius);
}

template<DirectionSampler DirectionSamplerType, IntervalSampler IntervalSamplerType>
bool HarmonicSpawnShell<DirectionSamplerType, IntervalSamplerType>::returns(
    double distance,
    double radius
) {
    return _interval_sampler.sample(UNIT_INTERVAL) < return_probability(distance, radius);
}

template<DirectionSampler DirectionSamplerType, IntervalSampler IntervalSamplerType>
Vector3 HarmonicSpawnShell<DirectionSamplerType, IntervalSamplerType>::landing(
    const Vector3& position,
    double distance,
    double radius
) {
    double cosine = harmonic_hit_cosine(distance, radius, _interval_sampler.sample(UNIT_INTERVAL));
    double azimuth = _interval_sampler.sample(Interval(0.0, TWO_PI));
    double ring_radius = std::sqrt(std::max(0.0, 1.0 - cosine * cosine));

    Vector3 axis = position / distance;
    Vector3 tangent = axis.unitOrthogonal();
    Vector3 bitangent = axis.cross(tangent);

    return radius * (
        cosine * axis +
        ring_radius * (std::cos(azimuth) * tangent + std::sin(azimuth) * bitangent)
    );
}

} // namespace geometry_art::dla

#endif //GEOMETRY_ART_DLA_SPAWN_SHELL_HARMONIC_SPAWN_SHELL_HPP_
