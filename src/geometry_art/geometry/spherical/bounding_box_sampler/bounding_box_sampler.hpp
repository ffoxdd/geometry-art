#ifndef GEOMETRY_ART_GEOMETRY_SPHERICAL_BOUNDING_BOX_SAMPLER_HPP_
#define GEOMETRY_ART_GEOMETRY_SPHERICAL_BOUNDING_BOX_SAMPLER_HPP_

#include "../../../types.hpp"
#include "../bounding_box.hpp"
#include <concepts>

namespace geometry_art::geometry::spherical {

using geometry_art::VectorS2;

template<typename T>
concept BoundingBoxSampler = requires(T sampler, const BoundingBox& bounding_box) {
    { sampler.sample(bounding_box) } -> std::same_as<VectorS2>;
};

} // namespace geometry_art::geometry::spherical

namespace geometry_art {
template<typename T>
concept SphericalBoundingBoxSampler = geometry::spherical::BoundingBoxSampler<T>;
}

#endif //GEOMETRY_ART_GEOMETRY_SPHERICAL_BOUNDING_BOX_SAMPLER_HPP_
