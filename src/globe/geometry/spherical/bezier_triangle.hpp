#ifndef GLOBEART_SRC_GLOBE_GEOMETRY_SPHERICAL_BEZIER_TRIANGLE_HPP_
#define GLOBEART_SRC_GLOBE_GEOMETRY_SPHERICAL_BEZIER_TRIANGLE_HPP_

#include "barycentric.hpp"
#include "../../types.hpp"
#include "../../math/polynomial/multi_index.hpp"
#include "../../math/polynomial/polynomial.hpp"
#include <CGAL/assertions.h>
#include <algorithm>
#include <array>
#include <cstddef>
#include <vector>

namespace globe::geometry::spherical {

using globe::math::polynomial::MultiIndex;
using globe::math::polynomial::Polynomial;

// A spherical polynomial written against the corners of a triangle rather
// than against the coordinate axes. Two properties are what it is for:
//
//   - the polynomial is bounded by its coefficients, so nonnegative
//     coefficients guarantee a nonnegative density and no sampling is needed
//     to know it;
//   - smoothness across a shared edge is a linear relation between the
//     coefficients of the two triangles, so a C1 field can be constructed
//     rather than fitted and hoped for.
//
// Each coefficient belongs to a domain point, the corner combination
// `(i, j, k) / degree`, which is where the interpolation conditions are
// written.
class BezierTriangle {
 public:
    BezierTriangle(const std::array<VectorS2, 3>& corners, int degree);

    [[nodiscard]] int degree() const { return _degree; }
    [[nodiscard]] const Barycentric& barycentric() const { return _barycentric; }

    [[nodiscard]] double coefficient(const MultiIndex& index) const;
    void set_coefficient(const MultiIndex& index, double value);

    [[nodiscard]] double lowest_coefficient() const;
    [[nodiscard]] Vector3 domain_point(const MultiIndex& index) const;
    [[nodiscard]] Polynomial to_polynomial() const;

    [[nodiscard]] static std::vector<MultiIndex> indices(int degree);
    [[nodiscard]] static double multinomial(const MultiIndex& index);

 private:
    Barycentric _barycentric;
    int _degree;
    std::vector<double> _coefficients;

    [[nodiscard]] size_t slot(const MultiIndex& index) const;
};

inline BezierTriangle::BezierTriangle(const std::array<VectorS2, 3>& corners, int degree) :
    _barycentric(corners),
    _degree(degree),
    _coefficients(indices(degree).size(), 0.0) {
    CGAL_precondition(degree >= 0);
}

inline double BezierTriangle::coefficient(const MultiIndex& index) const {
    return _coefficients[slot(index)];
}

inline void BezierTriangle::set_coefficient(const MultiIndex& index, double value) {
    _coefficients[slot(index)] = value;
}

inline double BezierTriangle::lowest_coefficient() const {
    return *std::min_element(_coefficients.begin(), _coefficients.end());
}

inline Vector3 BezierTriangle::domain_point(const MultiIndex& index) const {
    CGAL_precondition(index.degree() == _degree);

    return _barycentric.combination(
        Vector3(index.x, index.y, index.z) / static_cast<double>(_degree)
    );
}

// Each term is a product of the barycentric coordinates' linear forms, so
// the conversion is a sequence of linear-form multiplications rather than a
// change of basis solved numerically.
inline Polynomial BezierTriangle::to_polynomial() const {
    Polynomial result(_degree);

    for (const MultiIndex& index : indices(_degree)) {
        double coefficient = this->coefficient(index);

        if (coefficient == 0.0) {
            continue;
        }

        Polynomial term = Polynomial::constant(coefficient * multinomial(index));

        for (int axis = 0; axis < 3; ++axis) {
            for (int power = 0; power < index[axis]; ++power) {
                term = term.times_linear_form(_barycentric.linear_form(axis));
            }
        }

        result = result.plus(term);
    }

    return result;
}

inline std::vector<MultiIndex> BezierTriangle::indices(int degree) {
    return MultiIndex::all_of_degree(degree);
}

inline double BezierTriangle::multinomial(const MultiIndex& index) {
    auto factorial = [](int value) {
        double result = 1.0;

        for (int i = 2; i <= value; ++i) {
            result *= i;
        }

        return result;
    };

    return factorial(index.degree()) / (factorial(index.x) * factorial(index.y) * factorial(index.z));
}

// Mirrors the order `MultiIndex::all_of_degree` generates: first index
// descending, then the second.
inline size_t BezierTriangle::slot(const MultiIndex& index) const {
    CGAL_precondition(index.degree() == _degree);

    int remaining = _degree - index.x;
    return static_cast<size_t>(remaining * (remaining + 1) / 2 + (remaining - index.y));
}

} // namespace globe::geometry::spherical

#endif //GLOBEART_SRC_GLOBE_GEOMETRY_SPHERICAL_BEZIER_TRIANGLE_HPP_
