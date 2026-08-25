#ifndef GEOMETRY_ART_GEOMETRY_CARTESIAN_BOUNDING_BOX_SAMPLER_BOUNDING_BOX_SAMPLER_HPP_
#define GEOMETRY_ART_GEOMETRY_CARTESIAN_BOUNDING_BOX_SAMPLER_BOUNDING_BOX_SAMPLER_HPP_

#include "../../../types.hpp"
#include "../bounding_box.hpp"
#include <concepts>

namespace geometry_art::geometry::cartesian {

template<typename T>
concept BoundingBoxSampler = requires(T sampler, const BoundingBox& bounding_box) {
    { sampler.sample(bounding_box) } -> std::same_as<Vector3>;
};

} // namespace geometry_art::geometry::cartesian

namespace geometry_art {
template<typename T>
concept BoundingBoxSampler = geometry::cartesian::BoundingBoxSampler<T>;
}

#endif //GEOMETRY_ART_GEOMETRY_CARTESIAN_BOUNDING_BOX_SAMPLER_BOUNDING_BOX_SAMPLER_HPP_
