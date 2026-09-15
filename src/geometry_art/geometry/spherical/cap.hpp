#ifndef GEOMETRY_ART_GEOMETRY_SPHERICAL_CAP_HPP_
#define GEOMETRY_ART_GEOMETRY_SPHERICAL_CAP_HPP_

#include "../../types.hpp"
#include <algorithm>
#include <cmath>

namespace geometry_art::geometry::spherical {

// The points of the sphere on the near side of a plane: everything whose
// height along the axis is at least the offset. At offset zero it is a
// hemisphere, and moving the plane inward by a geodesic distance shrinks it
// to the cap whose rim runs parallel to the great circle that far inside.
struct Cap {
    VectorS2 axis;
    double offset;

    [[nodiscard]] static Cap hemisphere(const VectorS2& axis) { return Cap{axis, 0.0}; }
    [[nodiscard]] static Cap hemisphere_inset_by(const VectorS2& axis, double distance) { return Cap{axis, std::sin(distance)}; }

    [[nodiscard]] double excess(const VectorS2& point) const { return axis.dot(point) - offset; }
    [[nodiscard]] bool contains(const VectorS2& point) const { return excess(point) >= 0.0; }
    [[nodiscard]] VectorS2 rim_center() const { return offset * axis; }
    [[nodiscard]] double rim_radius() const { return std::sqrt(std::max(0.0, 1.0 - offset * offset)); }
};

} // namespace geometry_art::geometry::spherical

#endif //GEOMETRY_ART_GEOMETRY_SPHERICAL_CAP_HPP_
