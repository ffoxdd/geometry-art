#ifndef GLOBEART_SRC_GLOBE_GEOMETRY_SPHERICAL_BARYCENTRIC_HPP_
#define GLOBEART_SRC_GLOBE_GEOMETRY_SPHERICAL_BARYCENTRIC_HPP_

#include "../../types.hpp"
#include <CGAL/assertions.h>
#include <Eigen/Dense>
#include <array>

namespace globe::geometry::spherical {

// Coordinates of a point against the three corners of a spherical triangle,
// read as directions in space rather than as a partition of unity: a point
// on the sphere is the combination of the corners that reaches it, and the
// coordinates sum to one only at the corners themselves.
//
// They are linear in the point, which is what lets the whole Bernstein-Bezier
// apparatus -- the smoothness conditions above all -- carry over from the
// plane unchanged.
class Barycentric {
 public:
    explicit Barycentric(const std::array<VectorS2, 3>& corners);

    [[nodiscard]] Vector3 coordinates(const Vector3& point) const;
    [[nodiscard]] Vector3 combination(const Vector3& coordinates) const;
    [[nodiscard]] const std::array<VectorS2, 3>& corners() const { return _corners; }

    // The linear form reading off one coordinate, so that a coordinate can be
    // multiplied into a polynomial rather than only evaluated.
    [[nodiscard]] Vector3 linear_form(int index) const;

 private:
    std::array<VectorS2, 3> _corners;
    Matrix3 _to_coordinates;
    Matrix3 _to_point;
};

inline Barycentric::Barycentric(const std::array<VectorS2, 3>& corners) :
    _corners(corners) {

    _to_point.col(0) = corners[0];
    _to_point.col(1) = corners[1];
    _to_point.col(2) = corners[2];

    Eigen::FullPivLU<Matrix3> decomposition(_to_point);
    CGAL_precondition(decomposition.isInvertible());
    _to_coordinates = decomposition.inverse();
}

inline Vector3 Barycentric::coordinates(const Vector3& point) const {
    return _to_coordinates * point;
}

inline Vector3 Barycentric::combination(const Vector3& coordinates) const {
    return _to_point * coordinates;
}

inline Vector3 Barycentric::linear_form(int index) const {
    return _to_coordinates.row(index).transpose();
}

} // namespace globe::geometry::spherical

#endif //GLOBEART_SRC_GLOBE_GEOMETRY_SPHERICAL_BARYCENTRIC_HPP_
