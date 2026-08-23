#ifndef GLOBEART_SRC_GLOBE_FIELDS_SPHERICAL_LOCAL_QUADRATIC_FIT_HPP_
#define GLOBEART_SRC_GLOBE_FIELDS_SPHERICAL_LOCAL_QUADRATIC_FIT_HPP_

#include "../scalar/field.hpp"
#include "../../types.hpp"
#include "../../math/polynomial/multi_index.hpp"
#include <CGAL/assertions.h>
#include <Eigen/Dense>
#include <cmath>
#include <cstddef>
#include <vector>

namespace globe::fields::spherical {

using globe::math::polynomial::MultiIndex;

// Reads a field's value and gradient at a point by fitting a quadratic to
// many samples across a neighbourhood, rather than by evaluating the field
// and its difference quotient at the point itself.
//
// The difference is aliasing. A point sample carries whatever the field
// happens to do at that point, including structure far finer than the mesh
// that will represent it, and that structure then appears in the
// representation as though it were real at the mesh's scale -- corrupting
// the local average, which is the one thing a tessellation reads. A fit
// over a neighbourhood leaves sub-neighbourhood structure in the residual
// where it belongs. Reading a *gradient* at a point is the worst case of
// this: on a rough field a difference quotient describes the roughness and
// not the trend.
//
// The fit is onto the same space the pieces live in -- quadratics
// homogeneous in space -- so anything already in that space is returned
// exactly, and nothing is lost for fields that are smooth.
class LocalQuadraticFit {
 public:
    struct Reading {
        double value;
        Vector3 tangential_gradient;
    };

    LocalQuadraticFit(double radius, size_t sample_count = DEFAULT_SAMPLE_COUNT);

    template<scalar::Field ScalarFieldType>
    [[nodiscard]] Reading at(ScalarFieldType& scalar_field, const VectorS2& point) const;

    [[nodiscard]] std::vector<VectorS2> samples_around(const VectorS2& point) const;

 private:
    static constexpr size_t DEFAULT_SAMPLE_COUNT = 48;
    static constexpr size_t COEFFICIENT_COUNT = 6;

    double _radius;
    size_t _sample_count;

    [[nodiscard]] static Matrix3 quadratic_form(const Eigen::VectorXd& coefficients);
    [[nodiscard]] static std::vector<MultiIndex> basis();
};

inline LocalQuadraticFit::LocalQuadraticFit(double radius, size_t sample_count) :
    _radius(radius),
    _sample_count(sample_count) {
    CGAL_precondition(radius > 0.0);
    CGAL_precondition(sample_count >= COEFFICIENT_COUNT);
}

template<scalar::Field ScalarFieldType>
LocalQuadraticFit::Reading LocalQuadraticFit::at(ScalarFieldType& scalar_field, const VectorS2& point) const {
    std::vector<VectorS2> samples = samples_around(point);
    std::vector<MultiIndex> monomials = basis();

    Eigen::MatrixXd design(samples.size(), COEFFICIENT_COUNT);
    Eigen::VectorXd values(samples.size());

    for (size_t row = 0; row < samples.size(); ++row) {
        const VectorS2& sample = samples[row];
        values[static_cast<Eigen::Index>(row)] = scalar_field.value(sample);

        for (size_t column = 0; column < COEFFICIENT_COUNT; ++column) {
            const MultiIndex& index = monomials[column];
            design(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(column)) =
                std::pow(sample.x(), index.x) * std::pow(sample.y(), index.y) * std::pow(sample.z(), index.z);
        }
    }

    Matrix3 form = quadratic_form(design.colPivHouseholderQr().solve(values));
    Vector3 gradient = 2.0 * form * point;

    return Reading{
        point.dot(form * point),
        gradient - gradient.dot(point) * point
    };
}

// A spiral over the cap, equal-area in the polar angle so the samples do
// not crowd the centre and leave the rim to a handful of points.
inline std::vector<VectorS2> LocalQuadraticFit::samples_around(const VectorS2& point) const {
    Vector3 helper = std::abs(point.z()) < 0.9 ? Vector3::UnitZ() : Vector3::UnitX();
    Vector3 first = point.cross(helper).normalized();
    Vector3 second = point.cross(first);

    double lowest_cosine = std::cos(_radius);
    double golden_angle = M_PI * (3.0 - std::sqrt(5.0));

    std::vector<VectorS2> result;
    result.reserve(_sample_count);

    for (size_t index = 0; index < _sample_count; ++index) {
        double fraction = (static_cast<double>(index) + 0.5) / static_cast<double>(_sample_count);
        double cosine = 1.0 - fraction * (1.0 - lowest_cosine);
        double sine = std::sqrt(std::max(0.0, 1.0 - cosine * cosine));
        double azimuth = golden_angle * static_cast<double>(index);

        result.push_back(
            VectorS2(cosine * point + sine * (std::cos(azimuth) * first + std::sin(azimuth) * second)).normalized()
        );
    }

    return result;
}

inline Matrix3 LocalQuadraticFit::quadratic_form(const Eigen::VectorXd& coefficients) {
    std::vector<MultiIndex> monomials = basis();
    Matrix3 form = Matrix3::Zero();

    for (size_t column = 0; column < COEFFICIENT_COUNT; ++column) {
        const MultiIndex& index = monomials[column];
        double coefficient = coefficients[static_cast<Eigen::Index>(column)];

        int first = -1;
        int second = -1;

        for (int axis = 0; axis < 3; ++axis) {
            for (int repeat = 0; repeat < index[axis]; ++repeat) {
                (first < 0 ? first : second) = axis;
            }
        }

        if (first == second) {
            form(first, first) += coefficient;
        } else {
            form(first, second) += 0.5 * coefficient;
            form(second, first) += 0.5 * coefficient;
        }
    }

    return form;
}

inline std::vector<MultiIndex> LocalQuadraticFit::basis() {
    return MultiIndex::all_of_degree(2);
}

} // namespace globe::fields::spherical

#endif //GLOBEART_SRC_GLOBE_FIELDS_SPHERICAL_LOCAL_QUADRATIC_FIT_HPP_
