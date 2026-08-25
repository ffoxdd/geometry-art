#ifndef GEOMETRY_ART_GENERATORS_CARTESIAN_RANDOM_POINT_GENERATOR_HPP_
#define GEOMETRY_ART_GENERATORS_CARTESIAN_RANDOM_POINT_GENERATOR_HPP_

#include "../../types.hpp"
#include "../../geometry/cartesian/bounding_box.hpp"
#include "../../geometry/cartesian/bounding_box_sampler/bounding_box_sampler.hpp"
#include "../../geometry/cartesian/bounding_box_sampler/uniform_bounding_box_sampler.hpp"
#include <vector>

namespace geometry_art::generators::cartesian {

template<BoundingBoxSampler BoundingBoxSamplerType = UniformBoundingBoxSampler<>>
class RandomPointGenerator {
 public:
    RandomPointGenerator() = default;

    explicit RandomPointGenerator(BoundingBoxSamplerType sampler)
        : _sampler(std::move(sampler)) {
    }

    std::vector<Vector3> generate(size_t count);
    std::vector<Vector3> generate(size_t count, const BoundingBox &bounding_box);

 private:
    BoundingBoxSamplerType _sampler;
};

template<BoundingBoxSampler BoundingBoxSamplerType>
std::vector<Vector3>
RandomPointGenerator<BoundingBoxSamplerType>::generate(size_t count) {
    static const BoundingBox unit_cube = BoundingBox::unit_cube();
    return generate(count, unit_cube);
}

template<BoundingBoxSampler BoundingBoxSamplerType>
std::vector<Vector3>
RandomPointGenerator<BoundingBoxSamplerType>::generate(size_t count, const BoundingBox &bounding_box) {
    std::vector<Vector3> points;
    points.reserve(count);
    for (size_t i = 0; i < count; ++i) {
        points.push_back(_sampler.sample(bounding_box));
    }
    return points;
}

} // namespace geometry_art::generators::cartesian

#endif //GEOMETRY_ART_GENERATORS_CARTESIAN_RANDOM_POINT_GENERATOR_HPP_
