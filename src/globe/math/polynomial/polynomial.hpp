#ifndef GLOBEART_SRC_GLOBE_MATH_POLYNOMIAL_POLYNOMIAL_HPP_
#define GLOBEART_SRC_GLOBE_MATH_POLYNOMIAL_POLYNOMIAL_HPP_

#include "moments.hpp"
#include "monomial_table.hpp"
#include "multi_index.hpp"
#include "../../types.hpp"
#include <CGAL/assertions.h>
#include <algorithm>
#include <cmath>

namespace globe::math::polynomial {

class Polynomial {
 public:
    explicit Polynomial(int max_degree);

    [[nodiscard]] int max_degree() const { return _coefficients.max_degree(); }
    [[nodiscard]] const MonomialTable& coefficients() const { return _coefficients; }

    [[nodiscard]] double coefficient(const MultiIndex& index) const { return _coefficients.at(index); }
    void set_coefficient(const MultiIndex& index, double value) { _coefficients.set(index, value); }
    void add_coefficient(const MultiIndex& index, double value) { _coefficients.add(index, value); }

    [[nodiscard]] double value(const Vector3& point) const;
    [[nodiscard]] double integrate(const Moments& moments) const;
    [[nodiscard]] Polynomial times_coordinate(int axis) const;
    [[nodiscard]] Polynomial times_linear_form(const Vector3& form) const;
    [[nodiscard]] Polynomial plus(const Polynomial& other) const;
    [[nodiscard]] Polynomial scaled(double factor) const;

    [[nodiscard]] bool operator==(const Polynomial& other) const = default;

    [[nodiscard]] static Polynomial constant(double value);
    [[nodiscard]] static Polynomial linear(double constant, const Vector3& gradient);
    [[nodiscard]] static Polynomial quadratic(
        double constant,
        const Vector3& gradient,
        const Eigen::Matrix3d& quadratic_form
    );

 private:
    MonomialTable _coefficients;
};

inline Polynomial::Polynomial(int max_degree) :
    _coefficients(max_degree) {
}

inline double Polynomial::value(const Vector3& point) const {
    double result = 0.0;

    for (const MultiIndex& index : MultiIndex::all_up_to(max_degree())) {
        double coefficient = _coefficients.at(index);

        if (coefficient == 0.0) {
            continue;
        }

        result += coefficient *
            std::pow(point.x(), index.x) *
            std::pow(point.y(), index.y) *
            std::pow(point.z(), index.z);
    }

    return result;
}

inline double Polynomial::integrate(const Moments& moments) const {
    CGAL_precondition(moments.max_degree() >= max_degree());

    double result = 0.0;

    for (const MultiIndex& index : MultiIndex::all_up_to(max_degree())) {
        result += _coefficients.at(index) * moments.at(index);
    }

    return result;
}

inline Polynomial Polynomial::times_coordinate(int axis) const {
    Polynomial result(max_degree() + 1);

    for (const MultiIndex& index : MultiIndex::all_up_to(max_degree())) {
        result.set_coefficient(index.raised(axis), _coefficients.at(index));
    }

    return result;
}

// The generalisation of multiplying by one coordinate: any linear form
// `form . x`, which is what a barycentric coordinate is.
inline Polynomial Polynomial::times_linear_form(const Vector3& form) const {
    Polynomial result(max_degree() + 1);

    for (const MultiIndex& index : MultiIndex::all_up_to(max_degree())) {
        double coefficient = _coefficients.at(index);

        if (coefficient == 0.0) {
            continue;
        }

        for (int axis = 0; axis < 3; ++axis) {
            result.add_coefficient(index.raised(axis), coefficient * form[axis]);
        }
    }

    return result;
}

inline Polynomial Polynomial::plus(const Polynomial& other) const {
    Polynomial result(std::max(max_degree(), other.max_degree()));

    for (const MultiIndex& index : MultiIndex::all_up_to(max_degree())) {
        result.add_coefficient(index, _coefficients.at(index));
    }

    for (const MultiIndex& index : MultiIndex::all_up_to(other.max_degree())) {
        result.add_coefficient(index, other.coefficient(index));
    }

    return result;
}

inline Polynomial Polynomial::scaled(double factor) const {
    Polynomial result(max_degree());

    for (const MultiIndex& index : MultiIndex::all_up_to(max_degree())) {
        result.set_coefficient(index, _coefficients.at(index) * factor);
    }

    return result;
}

inline Polynomial Polynomial::constant(double value) {
    Polynomial result(0);
    result.set_coefficient(MultiIndex{0, 0, 0}, value);
    return result;
}

inline Polynomial Polynomial::linear(double constant, const Vector3& gradient) {
    Polynomial result(1);
    result.set_coefficient(MultiIndex{0, 0, 0}, constant);

    for (int axis = 0; axis < 3; ++axis) {
        result.set_coefficient(MultiIndex::unit(axis), gradient[axis]);
    }

    return result;
}

inline Polynomial Polynomial::quadratic(
    double constant,
    const Vector3& gradient,
    const Eigen::Matrix3d& quadratic_form
) {
    Polynomial result(2);
    result.set_coefficient(MultiIndex{0, 0, 0}, constant);

    for (int axis = 0; axis < 3; ++axis) {
        result.set_coefficient(MultiIndex::unit(axis), gradient[axis]);
    }

    for (int row = 0; row < 3; ++row) {
        for (int column = 0; column < 3; ++column) {
            result.add_coefficient(
                MultiIndex::unit(row).raised(column),
                quadratic_form(row, column)
            );
        }
    }

    return result;
}

} // namespace globe::math::polynomial

#endif //GLOBEART_SRC_GLOBE_MATH_POLYNOMIAL_POLYNOMIAL_HPP_
