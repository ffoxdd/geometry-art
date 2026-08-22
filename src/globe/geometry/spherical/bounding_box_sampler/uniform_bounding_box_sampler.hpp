#ifndef GLOBEART_SRC_GLOBE_GEOMETRY_SPHERICAL_UNIFORM_BOUNDING_BOX_SAMPLER_HPP_
#define GLOBEART_SRC_GLOBE_GEOMETRY_SPHERICAL_UNIFORM_BOUNDING_BOX_SAMPLER_HPP_

#include "../../../types.hpp"
#include "../bounding_box.hpp"
#include "../helpers.hpp"
#include "../../../math/interval_sampler/interval_sampler.hpp"
#include "../../../math/interval_sampler/uniform_interval_sampler.hpp"
#include "../../../math/circular_interval_sampler/circular_interval_sampler.hpp"
#include "../../../math/circular_interval_sampler/uniform_circular_interval_sampler.hpp"
#include <cmath>

namespace globe::geometry::spherical {

using globe::VectorS2;
using globe::IntervalSampler;
using globe::UniformIntervalSampler;
using globe::UniformCircularIntervalSampler;
using globe::CircularIntervalSampler;
using globe::TWO_PI;

template<
    IntervalSampler IntervalSamplerType = UniformIntervalSampler,
    CircularIntervalSampler<TWO_PI> CircularIntervalSamplerType = UniformCircularIntervalSampler
>
class UniformBoundingBoxSampler {
 public:
    UniformBoundingBoxSampler() = default;

    UniformBoundingBoxSampler(
        IntervalSamplerType interval_sampler,
        CircularIntervalSamplerType circular_interval_sampler
    ) : _interval_sampler(std::move(interval_sampler)),
        _circular_interval_sampler(std::move(circular_interval_sampler)) {}

    [[nodiscard]] inline VectorS2 sample(const BoundingBox &bounding_box) {
        double theta_val = _circular_interval_sampler.sample(bounding_box.theta_interval());
        double z = _interval_sampler.sample(bounding_box.z_interval());

        double r = std::sqrt(1.0 - z * z);
        return VectorS2(
            r * std::cos(theta_val),
            r * std::sin(theta_val),
            z
        );
    }

 private:
    IntervalSamplerType _interval_sampler;
    CircularIntervalSamplerType _circular_interval_sampler;
};

} // namespace globe::geometry::spherical

namespace globe {
template<
    IntervalSampler IntervalSamplerType = UniformIntervalSampler,
    CircularIntervalSampler<TWO_PI> CircularIntervalSamplerType = UniformCircularIntervalSampler
>
using UniformSphericalBoundingBoxSampler = geometry::spherical::UniformBoundingBoxSampler<IntervalSamplerType, CircularIntervalSamplerType>;
}

#endif //GLOBEART_SRC_GLOBE_GEOMETRY_SPHERICAL_UNIFORM_BOUNDING_BOX_SAMPLER_HPP_
