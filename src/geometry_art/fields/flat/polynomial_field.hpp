#ifndef GEOMETRY_ART_FIELDS_FLAT_POLYNOMIAL_FIELD_HPP_
#define GEOMETRY_ART_FIELDS_FLAT_POLYNOMIAL_FIELD_HPP_

#include "../region_integrals.hpp"
#include "../../types.hpp"
#include "../../geometry/planar/domain.hpp"
#include "../../geometry/planar/polygon.hpp"
#include "../../geometry/planar/segment.hpp"
#include "../../math/polynomial/moments.hpp"
#include "../../math/polynomial/multi_index.hpp"
#include "../../math/polynomial/polynomial.hpp"
#include <CGAL/assertions.h>
#include <array>
#include <utility>

namespace geometry_art::fields::flat {

using geometry::planar::Domain;
using geometry::planar::Polygon;
using geometry::planar::Segment;
using geometry_art::math::polynomial::Moments;
using geometry_art::math::polynomial::MultiIndex;
using geometry_art::math::polynomial::Polynomial;

// One polynomial density over the whole domain, integrated exactly through
// the regions' monomial moments. A polynomial is periodic along an axis only
// when it does not vary along it, so along every wrapped axis the density
// must be constant; a gradient needs a walled axis to rise along.
class PolynomialField {
 public:
    PolynomialField(Polynomial density, Domain domain);

    [[nodiscard]] static PolynomialField constant(double value, const Domain& domain);
    [[nodiscard]] static PolynomialField linear(double constant, const Vector2& gradient, const Domain& domain);

    [[nodiscard]] const Polynomial& density() const { return _density; }
    [[nodiscard]] const Domain& domain() const { return _domain; }
    [[nodiscard]] int degree() const { return _density.max_degree(); }
    [[nodiscard]] double total_mass() const { return _total_mass; }
    [[nodiscard]] double value(const Vector2& point) const;

    [[nodiscard]] RegionIntegrals integrals(const Polygon& polygon) const;
    [[nodiscard]] RegionIntegrals integrals(const Segment& segment) const;
    [[nodiscard]] double squared_norm_moment(const Polygon& polygon) const;
    [[nodiscard]] Matrix3 second_moment(const Segment& segment) const;
    [[nodiscard]] Vector3 gradient_masses(const Segment& segment) const;
    [[nodiscard]] Matrix3 gradient_first_moments(const Segment& segment) const;
    [[nodiscard]] std::array<Matrix3, 3> gradient_second_moments(const Segment& segment) const;

 private:
    Polynomial _density;
    Domain _domain;
    std::array<Polynomial, 2> _density_times_coordinate;
    Polynomial _density_times_squared_norm;
    std::array<Polynomial, 3> _density_times_coordinate_pair;
    std::array<Polynomial, 2> _derivatives;
    std::array<std::array<Polynomial, 2>, 2> _derivative_times_coordinate;
    std::array<std::array<Polynomial, 3>, 2> _derivative_times_coordinate_pair;
    double _total_mass;

    [[nodiscard]] bool constant_along_wrapped_axes() const;
    [[nodiscard]] RegionIntegrals integrals_from(const Moments& moments) const;
    [[nodiscard]] static Matrix3 planar_matrix(const std::array<Polynomial, 3>& pair_products, const Moments& moments);
    [[nodiscard]] static std::array<Polynomial, 2> coordinate_products(const Polynomial& polynomial);
    [[nodiscard]] static std::array<Polynomial, 3> pair_products(const std::array<Polynomial, 2>& products);
};

inline PolynomialField::PolynomialField(Polynomial density, Domain domain) :
    _density(std::move(density)),
    _domain(domain),
    _density_times_coordinate(coordinate_products(_density)),
    _density_times_squared_norm(
        _density_times_coordinate[0].times_coordinate(0).plus(_density_times_coordinate[1].times_coordinate(1))
    ),
    _density_times_coordinate_pair(pair_products(_density_times_coordinate)),
    _derivatives{_density.partial_derivative(0), _density.partial_derivative(1)},
    _derivative_times_coordinate{coordinate_products(_derivatives[0]), coordinate_products(_derivatives[1])},
    _derivative_times_coordinate_pair{
        pair_products(_derivative_times_coordinate[0]),
        pair_products(_derivative_times_coordinate[1])
    },
    _total_mass(0.0) {
    CGAL_precondition(constant_along_wrapped_axes());

    Polygon rectangle = Polygon::rectangle(Vector2::Zero(), Vector2(_domain.width, _domain.height));
    _total_mass = _density.integrate(rectangle.moments(degree()));
}

inline PolynomialField PolynomialField::constant(double value, const Domain& domain) {
    CGAL_precondition(value > 0.0);
    return PolynomialField(Polynomial::constant(value), domain);
}

inline PolynomialField PolynomialField::linear(double constant, const Vector2& gradient, const Domain& domain) {
    return PolynomialField(Polynomial::linear(constant, Vector3(gradient.x(), gradient.y(), 0.0)), domain);
}

inline double PolynomialField::value(const Vector2& point) const {
    return _density.value(Vector3(point.x(), point.y(), 0.0));
}

inline RegionIntegrals PolynomialField::integrals(const Polygon& polygon) const {
    return integrals_from(polygon.moments(degree() + 1));
}

inline RegionIntegrals PolynomialField::integrals(const Segment& segment) const {
    return integrals_from(segment.moments(degree() + 1));
}

inline double PolynomialField::squared_norm_moment(const Polygon& polygon) const {
    return _density_times_squared_norm.integrate(polygon.moments(degree() + 2));
}

inline Matrix3 PolynomialField::second_moment(const Segment& segment) const {
    return planar_matrix(_density_times_coordinate_pair, segment.moments(degree() + 2));
}

inline Vector3 PolynomialField::gradient_masses(const Segment& segment) const {
    Moments moments = segment.moments(degree());
    return Vector3(_derivatives[0].integrate(moments), _derivatives[1].integrate(moments), 0.0);
}

inline Matrix3 PolynomialField::gradient_first_moments(const Segment& segment) const {
    Moments moments = segment.moments(degree() + 1);
    Matrix3 result = Matrix3::Zero();

    for (int axis = 0; axis < 2; ++axis) {
        result(0, axis) = _derivative_times_coordinate[axis][0].integrate(moments);
        result(1, axis) = _derivative_times_coordinate[axis][1].integrate(moments);
    }

    return result;
}

inline std::array<Matrix3, 3> PolynomialField::gradient_second_moments(const Segment& segment) const {
    Moments moments = segment.moments(degree() + 2);

    return {
        planar_matrix(_derivative_times_coordinate_pair[0], moments),
        planar_matrix(_derivative_times_coordinate_pair[1], moments),
        Matrix3::Zero()
    };
}

inline bool PolynomialField::constant_along_wrapped_axes() const {
    for (int axis = 0; axis < 2; ++axis) {
        if (!_domain.wrapped(axis)) {
            continue;
        }

        for (const MultiIndex& index : MultiIndex::all_up_to(_derivatives[axis].max_degree())) {
            if (_derivatives[axis].coefficient(index) != 0.0) {
                return false;
            }
        }
    }

    return true;
}

inline RegionIntegrals PolynomialField::integrals_from(const Moments& moments) const {
    return RegionIntegrals{
        _density.integrate(moments),
        Vector3(
            _density_times_coordinate[0].integrate(moments),
            _density_times_coordinate[1].integrate(moments),
            0.0
        )
    };
}

inline Matrix3 PolynomialField::planar_matrix(const std::array<Polynomial, 3>& pair_products, const Moments& moments) {
    Matrix3 result = Matrix3::Zero();
    result(0, 0) = pair_products[0].integrate(moments);
    result(0, 1) = pair_products[1].integrate(moments);
    result(1, 0) = result(0, 1);
    result(1, 1) = pair_products[2].integrate(moments);
    return result;
}

inline std::array<Polynomial, 2> PolynomialField::coordinate_products(const Polynomial& polynomial) {
    return {polynomial.times_coordinate(0), polynomial.times_coordinate(1)};
}

inline std::array<Polynomial, 3> PolynomialField::pair_products(const std::array<Polynomial, 2>& products) {
    return {
        products[0].times_coordinate(0),
        products[0].times_coordinate(1),
        products[1].times_coordinate(1)
    };
}

} // namespace geometry_art::fields::flat

#endif //GEOMETRY_ART_FIELDS_FLAT_POLYNOMIAL_FIELD_HPP_
