#ifndef GLOBEART_SRC_GLOBE_FIELDS_SPHERICAL_HARMONIC_FIELD_HPP_
#define GLOBEART_SRC_GLOBE_FIELDS_SPHERICAL_HARMONIC_FIELD_HPP_

#include "field.hpp"
#include "../../types.hpp"
#include "../../geometry/spherical/arc.hpp"
#include "../../geometry/spherical/polygon/polygon.hpp"
#include <Eigen/Core>
#include <cmath>

namespace globe::fields::spherical {

using geometry::spherical::UNIT_SPHERE_AREA;

class HarmonicField {
 public:
    HarmonicField(
        double constant,
        const Eigen::Vector3d& linear,
        const Eigen::Matrix3d& quadratic
    );

    [[nodiscard]] double value(const VectorS2& point) const;
    [[nodiscard]] double mass(const Polygon& polygon) const;
    [[nodiscard]] double total_mass() const;
    [[nodiscard]] double edge_integral(const Arc& arc) const;
    [[nodiscard]] Eigen::Vector3d edge_gradient_integral(const Arc& arc) const;

 private:
    double _constant;
    Eigen::Vector3d _linear;
    Eigen::Matrix3d _quadratic;
    double _total_mass;

    [[nodiscard]] double evaluate_quadratic(const VectorS2& point) const;
    [[nodiscard]] double integrate_quadratic_polygon(const Polygon& polygon) const;
    [[nodiscard]] double integrate_quadratic_arc(const Arc& arc) const;
    [[nodiscard]] Eigen::Vector3d gradient_integrate_quadratic_arc(const Arc& arc) const;
};

inline HarmonicField::HarmonicField(
    double constant,
    const Eigen::Vector3d& linear,
    const Eigen::Matrix3d& quadratic
) :
    _constant(constant),
    _linear(linear),
    _quadratic(quadratic),
    _total_mass((constant + quadratic.trace() / 3.0) * UNIT_SPHERE_AREA) {
}

inline double HarmonicField::value(const VectorS2& point) const {
    return _constant + _linear.dot(point) + evaluate_quadratic(point);
}

inline double HarmonicField::evaluate_quadratic(const VectorS2& point) const {
    return point.transpose() * _quadratic * point;
}

inline double HarmonicField::mass(const Polygon& polygon) const {
    double constant_term = _constant * polygon.area();
    double linear_term = _linear.dot(polygon.first_moment());
    double quadratic_term = integrate_quadratic_polygon(polygon);

    return constant_term + linear_term + quadratic_term;
}

inline double HarmonicField::integrate_quadratic_polygon(const Polygon& polygon) const {
    Eigen::Matrix3d moment = polygon.second_moment();
    double result = 0.0;

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            result += _quadratic(i, j) * moment(i, j);
        }
    }

    return result;
}

inline double HarmonicField::total_mass() const {
    return _total_mass;
}

inline double HarmonicField::edge_integral(const Arc& arc) const {
    double constant_term = _constant * arc.length();
    double linear_term = _linear.dot(arc.first_moment());
    double quadratic_term = integrate_quadratic_arc(arc);

    return constant_term + linear_term + quadratic_term;
}

inline double HarmonicField::integrate_quadratic_arc(const Arc& arc) const {
    Eigen::Matrix3d moment = arc.second_moment();
    double result = 0.0;

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            result += _quadratic(i, j) * moment(i, j);
        }
    }

    return result;
}

inline Eigen::Vector3d HarmonicField::edge_gradient_integral(const Arc& arc) const {
    Eigen::Vector3d constant_term = _constant * arc.first_moment();
    Eigen::Vector3d linear_term = arc.second_moment() * _linear;
    Eigen::Vector3d quadratic_term = gradient_integrate_quadratic_arc(arc);

    return constant_term + linear_term + quadratic_term;
}

inline Eigen::Vector3d HarmonicField::gradient_integrate_quadratic_arc(const Arc& arc) const {
    constexpr size_t SAMPLES = 50;
    double arc_length = arc.length();

    if (arc_length < 1e-10) {
        return Eigen::Vector3d::Zero();
    }

    Eigen::Vector3d sum = Eigen::Vector3d::Zero();
    for (size_t i = 0; i < SAMPLES; ++i) {
        double t = (i + 0.5) / SAMPLES;
        VectorS2 point = arc.interpolate(t);
        double quad_value = evaluate_quadratic(point);
        sum += quad_value * point;
    }

    return sum * arc_length / SAMPLES;
}

static_assert(Field<HarmonicField>);

} // namespace globe::fields::spherical

#endif //GLOBEART_SRC_GLOBE_FIELDS_SPHERICAL_HARMONIC_FIELD_HPP_
