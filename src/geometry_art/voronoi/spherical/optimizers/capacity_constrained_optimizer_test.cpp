#include "capacity_constrained_optimizer.hpp"
#include "capacity_constrained_lagrangian.hpp"
#include "../core/random_builder.hpp"
#include "../../../fields/spherical/polynomial_field.hpp"
#include "../../../generators/spherical/fibonacci_point_generator.hpp"
#include "../../../testing/macros.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <memory>

using namespace geometry_art;
using namespace geometry_art::voronoi;
using namespace geometry_art::voronoi::spherical;
using fields::spherical::PolynomialField;
using generators::spherical::FibonacciPointGenerator;

namespace {

std::unique_ptr<Sphere> fibonacci_sphere(int count) {
    return RandomBuilder<FibonacciPointGenerator>(FibonacciPointGenerator()).build(count);
}

PolynomialField equator_dense_field() {
    Eigen::Matrix3d quadratic_form = Eigen::Matrix3d::Zero();
    quadratic_form(2, 2) = -0.9;
    return PolynomialField::quadratic(1.0, Vector3::Zero(), quadratic_form);
}

double relative_rms_capacity_error(const Sphere& sphere, const PolynomialField& field) {
    double target = field.total_mass() / sphere.size();
    double sum = 0.0;

    for (const auto& cell : sphere.cells()) {
        double error = field.integrals(cell).mass - target;
        sum += error * error;
    }

    return std::sqrt(sum / sphere.size()) / target;
}

double cvt_energy(const Sphere& sphere, const PolynomialField& field) {
    CapacityConstrainedLagrangian<> lagrangian(field, field.total_mass() / sphere.size());
    return lagrangian.evaluate(sphere, std::vector<double>(sphere.size(), 0.0), 0.0).cvt_energy;
}

}

TEST(CapacityConstrainedOptimizerTest, EqualizesCapacitiesForConstantField) {
    PolynomialField field = PolynomialField::constant(1.0);
    CapacityConstrainedParameters parameters;
    parameters.relative_capacity_tolerance = 1e-8;

    CapacityConstrainedOptimizer optimizer(fibonacci_sphere(8), field, parameters, noop_callback());
    auto sphere = optimizer.optimize();

    EXPECT_TRUE(optimizer.report().converged);
    EXPECT_LT(relative_rms_capacity_error(*sphere, field), 1e-8);
}

TEST(CapacityConstrainedOptimizerTest, EqualizesCapacitiesWithTheNewtonInnerSolver) {
    PolynomialField field = PolynomialField::constant(1.0);
    CapacityConstrainedParameters parameters;
    parameters.relative_capacity_tolerance = 1e-8;
    parameters.inner_solver = "newton";

    CapacityConstrainedOptimizer optimizer(fibonacci_sphere(8), field, parameters, noop_callback());
    auto sphere = optimizer.optimize();

    EXPECT_TRUE(optimizer.report().converged);
    EXPECT_LT(relative_rms_capacity_error(*sphere, field), 1e-8);
}

TEST(CapacityConstrainedOptimizerTest, EXPENSIVE_NewtonInnerSolverNeedsFewerIterationsThanLbfgs) {
    REQUIRE_EXPENSIVE();

    PolynomialField field = equator_dense_field();

    CapacityConstrainedParameters lbfgs_parameters;
    CapacityConstrainedOptimizer lbfgs(fibonacci_sphere(60), field, lbfgs_parameters, noop_callback());
    lbfgs.optimize();

    CapacityConstrainedParameters newton_parameters;
    newton_parameters.inner_solver = "newton";
    CapacityConstrainedOptimizer newton(fibonacci_sphere(60), field, newton_parameters, noop_callback());
    newton.optimize();

    EXPECT_TRUE(newton.report().converged);
    EXPECT_LT(newton.report().inner_iterations, lbfgs.report().inner_iterations);
}

TEST(CapacityConstrainedOptimizerTest, EqualizesCapacitiesForLinearField) {
    PolynomialField field = PolynomialField::linear(2.0, Vector3(0.0, 0.0, 1.0));
    CapacityConstrainedParameters parameters;
    parameters.relative_capacity_tolerance = 1e-8;

    CapacityConstrainedOptimizer optimizer(fibonacci_sphere(8), field, parameters, noop_callback());
    auto sphere = optimizer.optimize();

    EXPECT_TRUE(optimizer.report().converged);
    EXPECT_LT(relative_rms_capacity_error(*sphere, field), 1e-8);
}

TEST(CapacityConstrainedOptimizerTest, PreservesSiteCount) {
    PolynomialField field = PolynomialField::constant(1.0);
    CapacityConstrainedParameters parameters;
    parameters.max_outer_iterations = 1;
    parameters.max_inner_iterations = 5;

    CapacityConstrainedOptimizer optimizer(fibonacci_sphere(12), field, parameters, noop_callback());
    auto sphere = optimizer.optimize();

    EXPECT_EQ(sphere->size(), 12u);
}

TEST(CapacityConstrainedOptimizerTest, InvokesCallbackDuringOptimization) {
    PolynomialField field = PolynomialField::constant(1.0);
    CapacityConstrainedParameters parameters;
    parameters.max_outer_iterations = 1;
    parameters.max_inner_iterations = 3;
    size_t calls = 0;

    CapacityConstrainedOptimizer optimizer(fibonacci_sphere(8), field, parameters, [&](const Sphere&) { ++calls; });
    optimizer.optimize();

    EXPECT_GT(calls, 0u);
}

TEST(CapacityConstrainedOptimizerTest, EXPENSIVE_ConvergesForQuadraticFieldAndLowersCvtEnergy) {
    REQUIRE_EXPENSIVE();

    PolynomialField field = equator_dense_field();
    auto initial = fibonacci_sphere(60);
    double initial_energy = cvt_energy(*initial, field);
    CapacityConstrainedParameters parameters;
    parameters.relative_capacity_tolerance = 1e-7;

    CapacityConstrainedOptimizer optimizer(std::move(initial), field, parameters, noop_callback());
    auto sphere = optimizer.optimize();

    EXPECT_TRUE(optimizer.report().converged);
    EXPECT_LT(relative_rms_capacity_error(*sphere, field), 1e-7);
    EXPECT_LT(optimizer.report().cvt_energy, initial_energy);
}
