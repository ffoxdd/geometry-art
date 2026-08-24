#ifndef GLOBEART_SRC_GLOBE_FIELDS_FLAT_FIELD_HPP_
#define GLOBEART_SRC_GLOBE_FIELDS_FLAT_FIELD_HPP_

#include "../region_integrals.hpp"
#include "../../types.hpp"
#include <array>
#include "../../geometry/planar/polygon.hpp"
#include "../../geometry/planar/segment.hpp"
#include <concepts>

namespace globe::fields::flat {

using geometry::planar::Polygon;
using geometry::planar::Segment;

// A density on a flat domain. Every field is periodic in the domain's
// periods -- a cell reads the field in its own chart, and a value that
// depended on the chart would make the mass of a wrapped cell ambiguous.
template<typename T>
concept Field = requires(const T& field, const Vector2& point, const Polygon& polygon, const Segment& segment) {
    { field.value(point) } -> std::convertible_to<double>;
    { field.degree() } -> std::convertible_to<int>;
    { field.integrals(polygon) } -> std::convertible_to<RegionIntegrals>;
    { field.integrals(segment) } -> std::convertible_to<RegionIntegrals>;
    { field.squared_norm_moment(polygon) } -> std::convertible_to<double>;
    { field.second_moment(segment) } -> std::convertible_to<Matrix3>;
    { field.gradient_masses(segment) } -> std::convertible_to<Vector3>;
    { field.gradient_first_moments(segment) } -> std::convertible_to<Matrix3>;
    { field.gradient_second_moments(segment) } -> std::convertible_to<std::array<Matrix3, 3>>;
    { field.total_mass() } -> std::convertible_to<double>;
};

} // namespace globe::fields::flat

#endif //GLOBEART_SRC_GLOBE_FIELDS_FLAT_FIELD_HPP_
