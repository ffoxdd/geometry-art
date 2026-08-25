#ifndef GEOMETRY_ART_GEOMETRY_CARTESIAN_BOUNDING_BOX_SAMPLER_UNIFORM_BOUNDING_BOX_SAMPLER_HPP_
#define GEOMETRY_ART_GEOMETRY_CARTESIAN_BOUNDING_BOX_SAMPLER_UNIFORM_BOUNDING_BOX_SAMPLER_HPP_

#include "../../../types.hpp"
#include "../bounding_box.hpp"
#include "../../../math/interval_sampler/interval_sampler.hpp"
#include "../../../math/interval_sampler/uniform_interval_sampler.hpp"

namespace geometry_art::geometry::cartesian {

using geometry_art::IntervalSampler;
using geometry_art::UniformIntervalSampler;

template<IntervalSampler IntervalSamplerType = UniformIntervalSampler>
class UniformBoundingBoxSampler {
 public:
    UniformBoundingBoxSampler() = default;

    explicit UniformBoundingBoxSampler(IntervalSamplerType interval_sampler)
        : _interval_sampler(std::move(interval_sampler)) {}

    [[nodiscard]] inline Vector3 sample(const BoundingBox& bounding_box) {
        double x = _interval_sampler.sample(bounding_box.x_interval());
        double y = _interval_sampler.sample(bounding_box.y_interval());
        double z = _interval_sampler.sample(bounding_box.z_interval());

        return Vector3(x, y, z);
    }

 private:
    IntervalSamplerType _interval_sampler;
};

} // namespace geometry_art::geometry::cartesian

namespace geometry_art {
template<IntervalSampler IntervalSamplerType = UniformIntervalSampler>
using UniformBoundingBoxSampler = geometry::cartesian::UniformBoundingBoxSampler<IntervalSamplerType>;
}

#endif //GEOMETRY_ART_GEOMETRY_CARTESIAN_BOUNDING_BOX_SAMPLER_UNIFORM_BOUNDING_BOX_SAMPLER_HPP_
