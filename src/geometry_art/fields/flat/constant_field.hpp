#ifndef GEOMETRY_ART_FIELDS_FLAT_CONSTANT_FIELD_HPP_
#define GEOMETRY_ART_FIELDS_FLAT_CONSTANT_FIELD_HPP_

#include "../region_integrals.hpp"
#include "../../types.hpp"
#include "../../geometry/planar/polygon.hpp"
#include "../../geometry/planar/segment.hpp"
#include "../../math/polynomial/moments.hpp"
#include <CGAL/assertions.h>
#include <array>

namespace geometry_art::fields::flat {

using geometry::planar::Polygon;
using geometry::planar::Segment;
using geometry_art::math::polynomial::Moments;

// A uniform density on a rectangle of periods, the one field every flat
// domain admits without any question of periodicity.
class ConstantField {
 public:
    ConstantField(double density, double width, double height);

    [[nodiscard]] double value(const Vector2&) const { return _density; }
    [[nodiscard]] int degree() const { return 0; }
    [[nodiscard]] double total_mass() const { return _density * _width * _height; }

    [[nodiscard]] RegionIntegrals integrals(const Polygon& polygon) const;
    [[nodiscard]] RegionIntegrals integrals(const Segment& segment) const;
    [[nodiscard]] double squared_norm_moment(const Polygon& polygon) const;
    [[nodiscard]] Matrix3 second_moment(const Segment& segment) const;

    [[nodiscard]] Vector3 gradient_masses(const Segment&) const { return Vector3::Zero(); }
    [[nodiscard]] Matrix3 gradient_first_moments(const Segment&) const { return Matrix3::Zero(); }

    [[nodiscard]] std::array<Matrix3, 3> gradient_second_moments(const Segment&) const {
        return {Matrix3::Zero(), Matrix3::Zero(), Matrix3::Zero()};
    }

 private:
    double _density;
    double _width;
    double _height;

    [[nodiscard]] RegionIntegrals integrals_from(const Moments& moments) const;
    [[nodiscard]] Matrix3 second_moment_from(const Moments& moments) const;
};

inline ConstantField::ConstantField(double density, double width, double height) :
    _density(density),
    _width(width),
    _height(height) {
    CGAL_precondition(density > 0.0 && width > 0.0 && height > 0.0);
}

inline RegionIntegrals ConstantField::integrals(const Polygon& polygon) const {
    return integrals_from(polygon.moments(1));
}

inline RegionIntegrals ConstantField::integrals(const Segment& segment) const {
    return integrals_from(segment.moments(1));
}

inline double ConstantField::squared_norm_moment(const Polygon& polygon) const {
    Moments moments = polygon.moments(2);
    return _density * (moments.at(2, 0, 0) + moments.at(0, 2, 0));
}

inline Matrix3 ConstantField::second_moment(const Segment& segment) const {
    return second_moment_from(segment.moments(2));
}

inline RegionIntegrals ConstantField::integrals_from(const Moments& moments) const {
    return RegionIntegrals{
        _density * moments.at(0, 0, 0),
        _density * Vector3(moments.at(1, 0, 0), moments.at(0, 1, 0), 0.0)
    };
}

inline Matrix3 ConstantField::second_moment_from(const Moments& moments) const {
    Matrix3 result = Matrix3::Zero();
    result(0, 0) = moments.at(2, 0, 0);
    result(0, 1) = moments.at(1, 1, 0);
    result(1, 0) = moments.at(1, 1, 0);
    result(1, 1) = moments.at(0, 2, 0);
    return _density * result;
}

} // namespace geometry_art::fields::flat

#endif //GEOMETRY_ART_FIELDS_FLAT_CONSTANT_FIELD_HPP_
