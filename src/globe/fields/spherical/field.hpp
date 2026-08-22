#ifndef GLOBEART_SRC_GLOBE_FIELDS_SPHERICAL_FIELD_HPP_
#define GLOBEART_SRC_GLOBE_FIELDS_SPHERICAL_FIELD_HPP_

#include "region_integrals.hpp"
#include "../../types.hpp"
#include "../../geometry/spherical/arc.hpp"
#include "../../geometry/spherical/polygon/polygon.hpp"
#include "../../math/polynomial/moments.hpp"
#include <concepts>

namespace globe::fields::spherical {

using globe::math::polynomial::Moments;

template<typename T>
concept Field = requires(
    const T& field,
    const VectorS2& point,
    const Polygon& polygon,
    const Arc& arc,
    const Moments& moments
) {
    { field.value(point) } -> std::convertible_to<double>;
    { field.degree() } -> std::convertible_to<int>;
    { field.integrals(polygon) } -> std::convertible_to<RegionIntegrals>;
    { field.integrals(arc) } -> std::convertible_to<RegionIntegrals>;
    { field.integrals(moments) } -> std::convertible_to<RegionIntegrals>;
    { field.total_mass() } -> std::convertible_to<double>;
};

} // namespace globe::fields::spherical

#endif //GLOBEART_SRC_GLOBE_FIELDS_SPHERICAL_FIELD_HPP_
