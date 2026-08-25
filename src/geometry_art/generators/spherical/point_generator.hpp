#ifndef GEOMETRY_ART_GENERATORS_SPHERICAL_POINT_GENERATOR_HPP_
#define GEOMETRY_ART_GENERATORS_SPHERICAL_POINT_GENERATOR_HPP_

#include "../../types.hpp"
#include "../../geometry/spherical/bounding_box.hpp"
#include <vector>
#include <concepts>

namespace geometry_art::generators::spherical {

template<typename T>
concept PointGenerator = requires(T t, const SphericalBoundingBox &bounding_box, size_t count) {
    { t.generate(count) } -> std::convertible_to<std::vector<VectorS2>>;
    { t.generate(count, bounding_box) } -> std::convertible_to<std::vector<VectorS2>>;
    { t.last_attempt_count() } -> std::convertible_to<size_t>;
};

} // namespace geometry_art::generators::spherical

#endif //GEOMETRY_ART_GENERATORS_SPHERICAL_POINT_GENERATOR_HPP_
