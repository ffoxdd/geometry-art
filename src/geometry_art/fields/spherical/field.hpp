#ifndef GEOMETRY_ART_FIELDS_SPHERICAL_FIELD_HPP_
#define GEOMETRY_ART_FIELDS_SPHERICAL_FIELD_HPP_

#include "../region_integrals.hpp"
#include "../../types.hpp"
#include "../../geometry/spherical/arc.hpp"
#include "../../geometry/spherical/polygon/polygon.hpp"
#include "../../math/polynomial/moments.hpp"
#include <array>
#include <concepts>
#include <vector>

namespace geometry_art::fields::spherical {

using geometry_art::math::polynomial::Moments;

template<typename T>
concept Field = requires(
    const T& field,
    const VectorS2& point,
    const Polygon& polygon,
    const Arc& arc,
    const Moments& arc_moments,
    const std::vector<Moments>& polygon_arc_moments
) {
    { field.value(point) } -> std::convertible_to<double>;
    { field.degree() } -> std::convertible_to<int>;
    { field.integrals(polygon) } -> std::convertible_to<RegionIntegrals>;
    { field.integrals(arc) } -> std::convertible_to<RegionIntegrals>;
    { field.integrals(polygon, polygon_arc_moments) } -> std::convertible_to<RegionIntegrals>;
    { field.integrals(arc, arc_moments) } -> std::convertible_to<RegionIntegrals>;
    { field.second_moment(arc) } -> std::convertible_to<Matrix3>;
    { field.second_moment(arc, arc_moments) } -> std::convertible_to<Matrix3>;
    { field.gradient_second_moments(arc) } -> std::convertible_to<std::array<Matrix3, 3>>;
    { field.total_mass() } -> std::convertible_to<double>;
};

} // namespace geometry_art::fields::spherical

#endif //GEOMETRY_ART_FIELDS_SPHERICAL_FIELD_HPP_
