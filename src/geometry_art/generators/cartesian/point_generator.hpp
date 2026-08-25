#ifndef GEOMETRY_ART_GENERATORS_CARTESIAN_POINT_GENERATOR_HPP_
#define GEOMETRY_ART_GENERATORS_CARTESIAN_POINT_GENERATOR_HPP_

#include "../../types.hpp"
#include "../../geometry/cartesian/bounding_box.hpp"
#include <vector>
#include <concepts>

namespace geometry_art::generators::cartesian {

template<typename T>
concept PointGenerator = requires(T t, const BoundingBox &bounding_box, size_t count) {
    { t.generate(count) } -> std::convertible_to<std::vector<Vector3>>;
    { t.generate(count, bounding_box) } -> std::convertible_to<std::vector<Vector3>>;
};

} // namespace geometry_art::generators::cartesian

#endif //GEOMETRY_ART_GENERATORS_CARTESIAN_POINT_GENERATOR_HPP_
