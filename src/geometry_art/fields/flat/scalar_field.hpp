#ifndef GEOMETRY_ART_FIELDS_FLAT_SCALAR_FIELD_HPP_
#define GEOMETRY_ART_FIELDS_FLAT_SCALAR_FIELD_HPP_

#include "../../types.hpp"
#include <concepts>

namespace geometry_art::fields::flat {

// A pointwise density on the flat domain, the raw material a piecewise
// representation is sampled from. Implementations must be periodic in the
// domain's rectangle.
template<typename T>
concept ScalarField = requires(T field, const Vector2& point) {
    { field.value(point) } -> std::convertible_to<double>;
};

} // namespace geometry_art::fields::flat

#endif //GEOMETRY_ART_FIELDS_FLAT_SCALAR_FIELD_HPP_
