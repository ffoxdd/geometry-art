#ifndef GLOBEART_SRC_GLOBE_FIELDS_SPHERICAL_POLYNOMIAL_FIELD_FITTER_HPP_
#define GLOBEART_SRC_GLOBE_FIELDS_SPHERICAL_POLYNOMIAL_FIELD_FITTER_HPP_

#include "polynomial_field.hpp"
#include "../scalar/field.hpp"
#include "../../types.hpp"
#include "../../generators/spherical/point_generator.hpp"
#include "../../generators/spherical/fibonacci_point_generator.hpp"
#include "../../math/polynomial/multi_index.hpp"
#include "../../math/polynomial/polynomial.hpp"
#include <CGAL/assertions.h>
#include <Eigen/Dense>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

namespace globe::fields::spherical {

using globe::math::polynomial::MultiIndex;
using globe::math::polynomial::Polynomial;

template<generators::spherical::PointGenerator PointGeneratorType = generators::spherical::FibonacciPointGenerator>
class PolynomialFieldFitter {
 public:
    struct Fit {
        PolynomialField field;
        double root_mean_square_residual;

        // Least squares constrains no value, so a fit of a function with a
        // floor overshoots below it. A density that reaches zero cannot be
        // balanced against the others and one that goes negative is not a
        // density at all, so the lowest value the fit takes is part of the
        // result rather than something the caller has to go looking for.
        double lowest_sampled_value;
    };

    PolynomialFieldFitter(int degree, size_t sample_count, PointGeneratorType point_generator);

    template<scalar::Field ScalarFieldType>
    [[nodiscard]] Fit fit(ScalarFieldType& scalar_field);

 private:
    int _degree;
    size_t _sample_count;
    PointGeneratorType _point_generator;

    [[nodiscard]] std::vector<MultiIndex> basis() const;
    [[nodiscard]] static Eigen::MatrixXd design_matrix(
        const std::vector<VectorS2>& points,
        const std::vector<MultiIndex>& basis
    );
    [[nodiscard]] static double monomial(const VectorS2& point, const MultiIndex& index);
};

template<generators::spherical::PointGenerator PointGeneratorType>
PolynomialFieldFitter<PointGeneratorType>::PolynomialFieldFitter(
    int degree,
    size_t sample_count,
    PointGeneratorType point_generator
) :
    _degree(degree),
    _sample_count(sample_count),
    _point_generator(std::move(point_generator)) {
    CGAL_precondition(degree >= 0);
}

template<generators::spherical::PointGenerator PointGeneratorType>
template<scalar::Field ScalarFieldType>
typename PolynomialFieldFitter<PointGeneratorType>::Fit
PolynomialFieldFitter<PointGeneratorType>::fit(ScalarFieldType& scalar_field) {
    std::vector<VectorS2> points = _point_generator.generate(_sample_count);
    std::vector<MultiIndex> basis = this->basis();

    Eigen::VectorXd values(points.size());
    for (size_t i = 0; i < points.size(); ++i) {
        values[static_cast<Eigen::Index>(i)] = scalar_field.value(points[i]);
    }

    Eigen::MatrixXd design = design_matrix(points, basis);
    Eigen::VectorXd coefficients = design.colPivHouseholderQr().solve(values);
    double residual = (design * coefficients - values).norm() / std::sqrt(static_cast<double>(points.size()));

    Polynomial polynomial(_degree);
    for (size_t i = 0; i < basis.size(); ++i) {
        polynomial.set_coefficient(basis[i], coefficients[static_cast<Eigen::Index>(i)]);
    }

    Eigen::VectorXd fitted = design * coefficients;

    return Fit{PolynomialField(std::move(polynomial)), residual, fitted.minCoeff()};
}

template<generators::spherical::PointGenerator PointGeneratorType>
std::vector<MultiIndex> PolynomialFieldFitter<PointGeneratorType>::basis() const {
    std::vector<MultiIndex> result = MultiIndex::all_of_degree(_degree);

    if (_degree >= 1) {
        for (const MultiIndex& index : MultiIndex::all_of_degree(_degree - 1)) {
            result.push_back(index);
        }
    }

    return result;
}

template<generators::spherical::PointGenerator PointGeneratorType>
Eigen::MatrixXd PolynomialFieldFitter<PointGeneratorType>::design_matrix(
    const std::vector<VectorS2>& points,
    const std::vector<MultiIndex>& basis
) {
    Eigen::MatrixXd design(points.size(), basis.size());

    for (size_t row = 0; row < points.size(); ++row) {
        for (size_t column = 0; column < basis.size(); ++column) {
            design(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(column)) =
                monomial(points[row], basis[column]);
        }
    }

    return design;
}

template<generators::spherical::PointGenerator PointGeneratorType>
double PolynomialFieldFitter<PointGeneratorType>::monomial(const VectorS2& point, const MultiIndex& index) {
    return std::pow(point.x(), index.x) * std::pow(point.y(), index.y) * std::pow(point.z(), index.z);
}

} // namespace globe::fields::spherical

#endif //GLOBEART_SRC_GLOBE_FIELDS_SPHERICAL_POLYNOMIAL_FIELD_FITTER_HPP_
