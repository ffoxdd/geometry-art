#include "newton_optimizer.hpp"
#include "../lloyd_optimizer.hpp"
#include "../../../../fields/spherical/polynomial_field.hpp"
#include "../../../../testing/macros.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <memory>
#include <vector>

using namespace globe;
using namespace globe::voronoi::spherical;
using fields::spherical::PolynomialField;

namespace {

// Deliberately irregular so that both optimizers have real work to do.
std::vector<VectorS2> perturbed_sites(size_t count) {
    std::vector<VectorS2> sites;
    double golden_angle = M_PI * (std::sqrt(5.0) - 1.0);

    for (size_t i = 0; i < count; ++i) {
        double y = 1.0 - 2.0 * (i + 0.5) / count;
        double radius = std::sqrt(1.0 - y * y);
        double theta = golden_angle * i + 0.6 * std::sin(4.0 * i);
        double tilt = 0.25 * std::cos(3.0 * i);
        Vector3 point(std::cos(theta) * radius, y + tilt, std::sin(theta) * radius);
        sites.emplace_back(point.normalized());
    }

    return sites;
}

std::unique_ptr<Sphere> build_sphere(const std::vector<VectorS2>& sites) {
    auto sphere = std::make_unique<Sphere>();

    for (const VectorS2& site : sites) {
        sphere->insert(cgal::to_point(site));
    }

    return sphere;
}

PolynomialField test_field() {
    Matrix3 quadratic_form = Matrix3::Zero();
    quadratic_form(2, 2) = -0.6;
    return PolynomialField::quadratic(1.0, Vector3(0.2, -0.1, 0.3), quadratic_form);
}

double cvt_energy(const Sphere& sphere, const PolynomialField& field) {
    CapacityConstrainedLagrangian<PolynomialField> lagrangian(field, field.total_mass() / sphere.size());
    return lagrangian.evaluate(sphere, std::vector<double>(sphere.size(), 0.0), 0.0).cvt_energy;
}

} // namespace

TEST(NewtonOptimizerTest, LowersTheEnergyOnAConstantField) {
    PolynomialField field = PolynomialField::constant(1.0);
    auto sphere = build_sphere(perturbed_sites(12));
    double before = cvt_energy(*sphere, field);

    NewtonOptimizer<PolynomialField> optimizer(std::move(sphere), field, NewtonParameters(), noop_callback());
    auto optimized = optimizer.optimize();

    EXPECT_LT(optimizer.report().cvt_energy, before);
    EXPECT_GT(optimizer.report().accepted_steps, 0u);
    EXPECT_EQ(optimized->size(), 12u);
}

TEST(NewtonOptimizerTest, DrivesTheGradientToZero) {
    PolynomialField field = test_field();
    auto sphere = build_sphere(perturbed_sites(16));

    NewtonOptimizer<PolynomialField> optimizer(std::move(sphere), field, NewtonParameters(), noop_callback());
    optimizer.optimize();

    EXPECT_TRUE(optimizer.report().converged);
    EXPECT_LT(optimizer.report().gradient_norm, 1e-7);
}

TEST(NewtonOptimizerTest, EXPENSIVE_ReachesLowerEnergyThanLloydInFewerPasses) {
    REQUIRE_EXPENSIVE();

    PolynomialField field = test_field();
    size_t sites = 60;

    LloydOptimizer<PolynomialField> lloyd(build_sphere(perturbed_sites(sites)), field, 50, noop_callback());
    auto relaxed = lloyd.optimize();
    double lloyd_energy = cvt_energy(*relaxed, field);

    NewtonParameters parameters;
    parameters.max_iterations = 50;
    NewtonOptimizer<PolynomialField> newton(
        build_sphere(perturbed_sites(sites)),
        field,
        parameters,
        noop_callback()
    );
    newton.optimize();

    EXPECT_LT(newton.report().cvt_energy, lloyd_energy);
    EXPECT_LT(newton.report().iterations, 50u);
    EXPECT_TRUE(newton.report().converged);
}
