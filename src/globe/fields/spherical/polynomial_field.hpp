#ifndef GLOBEART_SRC_GLOBE_FIELDS_SPHERICAL_POLYNOMIAL_FIELD_HPP_
#define GLOBEART_SRC_GLOBE_FIELDS_SPHERICAL_POLYNOMIAL_FIELD_HPP_

#include "field.hpp"
#include "region_integrals.hpp"
#include "../../types.hpp"
#include "../../geometry/spherical/arc.hpp"
#include "../../geometry/spherical/polygon/polygon.hpp"
#include "../../math/polynomial/moments.hpp"
#include "../../math/polynomial/polynomial.hpp"
#include <CGAL/assertions.h>
#include <array>
#include <utility>
#include <vector>

namespace globe::fields::spherical {

using globe::math::polynomial::Moments;
using globe::math::polynomial::Polynomial;

class PolynomialField {
 public:
    explicit PolynomialField(Polynomial density);

    [[nodiscard]] const Polynomial& density() const { return _density; }
    [[nodiscard]] int degree() const { return _density.max_degree(); }

    [[nodiscard]] double value(const VectorS2& point) const;
    [[nodiscard]] RegionIntegrals integrals(const Polygon& polygon) const;
    [[nodiscard]] RegionIntegrals integrals(const Arc& arc) const;
    [[nodiscard]] RegionIntegrals integrals(const Polygon& polygon, const std::vector<Moments>& arc_moments) const;
    [[nodiscard]] RegionIntegrals integrals(const Arc& arc, const Moments& arc_moments) const;
    [[nodiscard]] Matrix3 second_moment(const Arc& arc) const;
    [[nodiscard]] Matrix3 second_moment(const Arc& arc, const Moments& arc_moments) const;
    [[nodiscard]] static int second_moment_degree(int density_degree) { return density_degree + 2; }
    [[nodiscard]] double total_mass() const { return _total_mass; }

    [[nodiscard]] static PolynomialField constant(double value);
    [[nodiscard]] static PolynomialField linear(double constant, const Vector3& gradient);
    [[nodiscard]] static PolynomialField quadratic(
        double constant,
        const Vector3& gradient,
        const Eigen::Matrix3d& quadratic_form
    );

 private:
    Polynomial _density;
    std::array<Polynomial, 3> _density_times_coordinate;
    std::vector<Polynomial> _density_times_coordinate_pair;
    double _total_mass;

    [[nodiscard]] RegionIntegrals integrate(const Moments& moments) const;
    [[nodiscard]] Matrix3 second_moment_from(const Moments& moments) const;
    [[nodiscard]] static std::vector<Polynomial> coordinate_pair_products(const Polynomial& density);
};

inline PolynomialField::PolynomialField(Polynomial density) :
    _density(std::move(density)),
    _density_times_coordinate{
        _density.times_coordinate(0),
        _density.times_coordinate(1),
        _density.times_coordinate(2)
    },
    _density_times_coordinate_pair(coordinate_pair_products(_density)),
    _total_mass(_density.integrate(Moments::unit_sphere(_density.max_degree()))) {
}

inline double PolynomialField::value(const VectorS2& point) const {
    return _density.value(point);
}

inline RegionIntegrals PolynomialField::integrals(const Polygon& polygon) const {
    return integrate(polygon.moments(degree() + 1));
}

inline RegionIntegrals PolynomialField::integrals(const Arc& arc) const {
    return integrate(arc.moments(degree() + 1));
}

inline RegionIntegrals PolynomialField::integrals(const Polygon& polygon, const std::vector<Moments>& arc_moments) const {
    return integrate(polygon.moments(degree() + 1, arc_moments));
}

inline RegionIntegrals PolynomialField::integrals([[maybe_unused]] const Arc& arc, const Moments& arc_moments) const {
    return integrate(arc_moments);
}

inline Matrix3 PolynomialField::second_moment(const Arc& arc) const {
    return second_moment_from(arc.moments(second_moment_degree(degree())));
}

inline Matrix3 PolynomialField::second_moment(
    [[maybe_unused]] const Arc& arc,
    const Moments& arc_moments
) const {
    CGAL_precondition(arc_moments.max_degree() >= second_moment_degree(degree()));
    return second_moment_from(arc_moments);
}

inline Matrix3 PolynomialField::second_moment_from(const Moments& moments) const {
    Matrix3 result;

    for (int row = 0; row < 3; ++row) {
        for (int column = 0; column < 3; ++column) {
            result(row, column) = _density_times_coordinate_pair[row * 3 + column].integrate(moments);
        }
    }

    return result;
}

inline std::vector<Polynomial> PolynomialField::coordinate_pair_products(const Polynomial& density) {
    std::vector<Polynomial> products;
    products.reserve(9);

    for (int row = 0; row < 3; ++row) {
        Polynomial row_product = density.times_coordinate(row);

        for (int column = 0; column < 3; ++column) {
            products.push_back(row_product.times_coordinate(column));
        }
    }

    return products;
}

inline RegionIntegrals PolynomialField::integrate(const Moments& moments) const {
    return RegionIntegrals{
        _density.integrate(moments),
        Vector3(
            _density_times_coordinate[0].integrate(moments),
            _density_times_coordinate[1].integrate(moments),
            _density_times_coordinate[2].integrate(moments)
        )
    };
}

inline PolynomialField PolynomialField::constant(double value) {
    return PolynomialField(Polynomial::constant(value));
}

inline PolynomialField PolynomialField::linear(double constant, const Vector3& gradient) {
    return PolynomialField(Polynomial::linear(constant, gradient));
}

inline PolynomialField PolynomialField::quadratic(
    double constant,
    const Vector3& gradient,
    const Eigen::Matrix3d& quadratic_form
) {
    return PolynomialField(Polynomial::quadratic(constant, gradient, quadratic_form));
}

static_assert(Field<PolynomialField>);

} // namespace globe::fields::spherical

#endif //GLOBEART_SRC_GLOBE_FIELDS_SPHERICAL_POLYNOMIAL_FIELD_HPP_
