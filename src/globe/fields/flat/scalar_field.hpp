#ifndef GLOBEART_SRC_GLOBE_FIELDS_FLAT_SCALAR_FIELD_HPP_
#define GLOBEART_SRC_GLOBE_FIELDS_FLAT_SCALAR_FIELD_HPP_

#include "../../types.hpp"
#include <concepts>

namespace globe::fields::flat {

// A pointwise density on the flat domain, the raw material a piecewise
// representation is sampled from. Implementations must be periodic in the
// domain's rectangle.
template<typename T>
concept ScalarField = requires(T field, const Vector2& point) {
    { field.value(point) } -> std::convertible_to<double>;
};

} // namespace globe::fields::flat

#endif //GLOBEART_SRC_GLOBE_FIELDS_FLAT_SCALAR_FIELD_HPP_
