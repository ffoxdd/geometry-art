#include "geometry_art/fields/scalar/noise_field.hpp"
#include "geometry_art/fields/spherical/piecewise_polynomial_field.hpp"
#include "geometry_art/fields/spherical/polynomial_field_fitter.hpp"
#include "geometry_art/generators/spherical/fibonacci_point_generator.hpp"
#include "geometry_art/geometry/spherical/triangle_mesh.hpp"
#include "geometry_art/math/interval.hpp"
#include "geometry_art/testing/macros.hpp"
#include <gtest/gtest.h>
#include <algorithm>
#include <iostream>
#include <vector>

using namespace geometry_art;
using fields::scalar::NoiseField;
using fields::spherical::PiecewisePolynomialField;
using fields::spherical::PolynomialFieldFitter;
using geometry::spherical::TriangleMesh;

namespace {

constexpr double DENSITY_FLOOR = 0.2;
constexpr int FIT_DEGREE = 8;
constexpr size_t FIT_SAMPLES = 20000;
constexpr int MESH_SUBDIVISIONS = 4;
constexpr int MESH_DEGREE = 2;

struct Extremes {
    double lowest;
    double highest;
    VectorS2 lowest_at;
};

template<typename FieldType>
Extremes extremes_of(const FieldType& field, size_t sample_count) {
    Extremes result{std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(), VectorS2()};

    for (const VectorS2& point : generators::spherical::FibonacciPointGenerator().generate(sample_count)) {
        double value = field.value(point);

        if (value < result.lowest) {
            result.lowest = value;
            result.lowest_at = point;
        }

        result.highest = std::max(result.highest, value);
    }

    return result;
}

} // namespace

// The density the solver sees must stay positive everywhere: a cell whose
// region has no mass cannot be balanced against the others, so the optimizer
// grows it without bound.
TEST(FieldPositivityTest, EXPENSIVE_PiecewiseNoiseStaysAboveTheFloor) {
    REQUIRE_EXPENSIVE();
    NoiseField noise(Interval(DENSITY_FLOOR, 1.0));
    PiecewisePolynomialField field = PiecewisePolynomialField::sample(
        TriangleMesh::icosphere(MESH_SUBDIVISIONS), MESH_DEGREE, noise);

    Extremes extremes = extremes_of(field, 200000);
    std::cout << "piecewise: lowest " << extremes.lowest << " highest " << extremes.highest << std::endl;

    EXPECT_GT(extremes.lowest, 0.0);
}

// Least squares constrains no value, so fitting a floored function
// overshoots below the floor: the global path produces something that is
// not a density at all. This is why the piecewise representation is the one
// the solver uses, and the fitter reports the shortfall rather than hiding
// it.
TEST(FieldPositivityTest, EXPENSIVE_GlobalFitOfFlooredNoiseGoesNegative) {
    REQUIRE_EXPENSIVE();
    NoiseField noise(Interval(DENSITY_FLOOR, 1.0));
    PolynomialFieldFitter<> fitter(FIT_DEGREE, FIT_SAMPLES, generators::spherical::FibonacciPointGenerator());
    auto fit = fitter.fit(noise);

    Extremes extremes = extremes_of(fit.field, 200000);

    EXPECT_LT(extremes.lowest, 0.0);
    EXPECT_LE(fit.lowest_sampled_value, 0.0);
    EXPECT_NEAR(fit.lowest_sampled_value, extremes.lowest, 0.1);
}
