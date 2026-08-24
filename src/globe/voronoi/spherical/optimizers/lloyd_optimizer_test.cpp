#include "lloyd_optimizer.hpp"
#include "../core/random_builder.hpp"
#include "../../../fields/spherical/polynomial_field.hpp"
#include "../../../generators/spherical/fibonacci_point_generator.hpp"
#include <gtest/gtest.h>
#include <cmath>

using namespace globe;
using namespace globe::voronoi;
using namespace globe::voronoi::spherical;
using fields::spherical::PolynomialField;
using generators::spherical::FibonacciPointGenerator;

namespace {

double weighted_deviation(const Sphere& sphere, const PolynomialField& field) {
    double total = 0.0;
    size_t index = 0;

    for (const auto& cell : sphere.cells()) {
        VectorS2 centroid = field.integrals(cell).first_moment.normalized();
        double deviation = distance(to_vector_s2(sphere.site(index)), centroid);
        total += deviation * deviation;
        ++index;
    }

    return std::sqrt(total / sphere.size());
}

}

TEST(LloydOptimizerTest, ReducesWeightedDeviation) {
    PolynomialField field = PolynomialField::linear(2.0, Vector3(0.0, 0.0, 1.0));
    auto sphere = RandomBuilder().build(10);
    double initial = weighted_deviation(*sphere, field);

    LloydOptimizer optimizer(std::move(sphere), field, 5, noop_callback());
    auto optimized = optimizer.optimize();

    EXPECT_LT(optimizer.final_deviation(), initial);
    EXPECT_NEAR(optimizer.final_deviation(), weighted_deviation(*optimized, field), 1e-12);
}

TEST(LloydOptimizerTest, PreservesPointCount) {
    auto sphere = RandomBuilder().build(15);

    LloydOptimizer optimizer(std::move(sphere), PolynomialField::constant(1.0), 3, noop_callback());
    auto optimized = optimizer.optimize();

    EXPECT_EQ(optimized->size(), 15u);
}

TEST(LloydOptimizerTest, ZeroPassesLeavesSitesUntouched) {
    auto sphere = RandomBuilder<FibonacciPointGenerator>(FibonacciPointGenerator()).build(8);
    VectorS2 first_site = to_vector_s2(sphere->site(0));

    LloydOptimizer optimizer(std::move(sphere), PolynomialField::constant(1.0), 0, noop_callback());
    auto optimized = optimizer.optimize();

    EXPECT_NEAR((to_vector_s2(optimized->site(0)) - first_site).norm(), 0.0, 1e-15);
}
