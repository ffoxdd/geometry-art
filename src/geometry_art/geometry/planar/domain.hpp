#ifndef GEOMETRY_ART_GEOMETRY_PLANAR_DOMAIN_HPP_
#define GEOMETRY_ART_GEOMETRY_PLANAR_DOMAIN_HPP_

#include "../../types.hpp"
#include <CGAL/assertions.h>
#include <cmath>

namespace geometry_art::geometry::planar {

// How a rectangle's axis closes: wrapped identifies its two ends, walled
// bounds the domain there.
enum class Closure {
    wrapped,
    walled
};

// A rectangle whose two axes each close in their own way: the flat torus
// wraps both, the plane walls both, and the cylinder wraps its width and
// walls its height. Everything in the flat family reads its periodicity
// from here.
struct Domain {
    double width;
    double height;
    Closure across;
    Closure along;

    [[nodiscard]] static Domain torus(double width, double height);
    [[nodiscard]] static Domain cylinder(double width, double height);
    [[nodiscard]] static Domain plane(double width, double height);

    [[nodiscard]] double extent(int axis) const { return axis == 0 ? width : height; }
    [[nodiscard]] Closure closure(int axis) const { return axis == 0 ? across : along; }
    [[nodiscard]] bool wrapped(int axis) const { return closure(axis) == Closure::wrapped; }
    [[nodiscard]] bool has_wrapped_axis() const { return wrapped(0) || wrapped(1); }
    [[nodiscard]] bool has_walled_axis() const { return !wrapped(0) || !wrapped(1); }
    [[nodiscard]] double area() const { return width * height; }

    // The representative of a point's class: wrapped coordinates reduced
    // into the rectangle, walled ones left where they are.
    [[nodiscard]] Vector2 canonical(const Vector2& point) const;
};

inline Domain Domain::torus(double width, double height) {
    CGAL_precondition(width > 0.0 && height > 0.0);
    return Domain{width, height, Closure::wrapped, Closure::wrapped};
}

inline Domain Domain::cylinder(double width, double height) {
    CGAL_precondition(width > 0.0 && height > 0.0);
    return Domain{width, height, Closure::wrapped, Closure::walled};
}

inline Domain Domain::plane(double width, double height) {
    CGAL_precondition(width > 0.0 && height > 0.0);
    return Domain{width, height, Closure::walled, Closure::walled};
}

inline Vector2 Domain::canonical(const Vector2& point) const {
    Vector2 result = point;

    for (int axis = 0; axis < 2; ++axis) {
        if (!wrapped(axis)) {
            continue;
        }

        double reduced = std::fmod(point[axis], extent(axis));

        if (reduced < 0.0) {
            reduced += extent(axis);
        }

        result[axis] = reduced;
    }

    return result;
}

} // namespace geometry_art::geometry::planar

#endif //GEOMETRY_ART_GEOMETRY_PLANAR_DOMAIN_HPP_
