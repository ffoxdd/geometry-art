#ifndef GEOMETRY_ART_DLA_DIRECTION_SAMPLER_UNIFORM_DIRECTION_SAMPLER_HPP_
#define GEOMETRY_ART_DLA_DIRECTION_SAMPLER_UNIFORM_DIRECTION_SAMPLER_HPP_

#include "../../types.hpp"
#include "../../math/interval.hpp"
#include "../../math/interval_sampler/interval_sampler.hpp"
#include "../../math/interval_sampler/uniform_interval_sampler.hpp"
#include <cmath>
#include <utility>

namespace geometry_art::dla {

// Uniform on the unit sphere by Archimedes: the height is uniform, the
// azimuth independent.
template<IntervalSampler IntervalSamplerType = UniformIntervalSampler>
class UniformDirectionSampler {
 public:
    UniformDirectionSampler() = default;

    explicit UniformDirectionSampler(IntervalSamplerType interval_sampler) :
        _interval_sampler(std::move(interval_sampler)) {
    }

    [[nodiscard]] Vector3 sample();

 private:
    IntervalSamplerType _interval_sampler;
};

template<IntervalSampler IntervalSamplerType>
Vector3 UniformDirectionSampler<IntervalSamplerType>::sample() {
    double height = _interval_sampler.sample(Interval(-1.0, 1.0));
    double azimuth = _interval_sampler.sample(Interval(0.0, TWO_PI));
    double ring_radius = std::sqrt(1.0 - height * height);

    return Vector3(
        ring_radius * std::cos(azimuth),
        ring_radius * std::sin(azimuth),
        height
    );
}

} // namespace geometry_art::dla

#endif //GEOMETRY_ART_DLA_DIRECTION_SAMPLER_UNIFORM_DIRECTION_SAMPLER_HPP_
