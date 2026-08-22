#ifndef GLOBEART_SRC_GLOBE_FIELDS_SPHERICAL_POLYNOMIAL_FIELD_HPP_
#define GLOBEART_SRC_GLOBE_FIELDS_SPHERICAL_POLYNOMIAL_FIELD_HPP_

#include "field.hpp"
#include "region_integrals.hpp"
#include "../../types.hpp"
#include "../../geometry/spherical/arc.hpp"
#include "../../geometry/spherical/polygon/polygon.hpp"
#include "../../math/polynomial/moments.hpp"
#include "../../math/polynomial/polynomial.hpp"
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
    double _total_mass;

    [[nodiscard]] RegionIntegrals integrate(const Moments& moments) const;
};

inline PolynomialField::PolynomialField(Polynomial density) :
    _density(std::move(density)),
    _density_times_coordinate{
        _density.times_coordinate(0),
        _density.times_coordinate(1),
        _density.times_coordinate(2)
    },
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
