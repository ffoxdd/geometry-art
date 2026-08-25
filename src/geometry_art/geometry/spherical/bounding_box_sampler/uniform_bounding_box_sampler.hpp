#ifndef GEOMETRY_ART_GEOMETRY_SPHERICAL_UNIFORM_BOUNDING_BOX_SAMPLER_HPP_
#define GEOMETRY_ART_GEOMETRY_SPHERICAL_UNIFORM_BOUNDING_BOX_SAMPLER_HPP_

#include "../../../types.hpp"
#include "../bounding_box.hpp"
#include "../helpers.hpp"
#include "../../../math/interval_sampler/interval_sampler.hpp"
#include "../../../math/interval_sampler/uniform_interval_sampler.hpp"
#include "../../../math/circular_interval_sampler/circular_interval_sampler.hpp"
#include "../../../math/circular_interval_sampler/uniform_circular_interval_sampler.hpp"
#include <cmath>

namespace geometry_art::geometry::spherical {

using geometry_art::VectorS2;
using geometry_art::IntervalSampler;
using geometry_art::UniformIntervalSampler;
using geometry_art::UniformCircularIntervalSampler;
using geometry_art::CircularIntervalSampler;
using geometry_art::TWO_PI;

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

} // namespace geometry_art::geometry::spherical

namespace geometry_art {
template<
    IntervalSampler IntervalSamplerType = UniformIntervalSampler,
    CircularIntervalSampler<TWO_PI> CircularIntervalSamplerType = UniformCircularIntervalSampler
>
using UniformSphericalBoundingBoxSampler = geometry::spherical::UniformBoundingBoxSampler<IntervalSamplerType, CircularIntervalSamplerType>;
}

#endif //GEOMETRY_ART_GEOMETRY_SPHERICAL_UNIFORM_BOUNDING_BOX_SAMPLER_HPP_
