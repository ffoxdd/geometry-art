#include "capacity_constrained_lagrangian.hpp"
#include "../../../fields/spherical/polynomial_field.hpp"
#include "../../../testing/macros.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <memory>
#include <vector>

using namespace globe;
using namespace globe::voronoi::spherical;
using fields::spherical::PolynomialField;

namespace {

std::vector<VectorS2> fibonacci_sites(size_t count) {
    std::vector<VectorS2> sites;
    double golden_angle = M_PI * (std::sqrt(5.0) - 1.0);

    for (size_t i = 0; i < count; ++i) {
        double y = 1.0 - 2.0 * (i + 0.5) / count;
        double radius = std::sqrt(1.0 - y * y);
        double theta = golden_angle * i;
        sites.emplace_back(std::cos(theta) * radius, y, std::sin(theta) * radius);
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
    Eigen::Matrix3d quadratic_form = Eigen::Matrix3d::Zero();
    quadratic_form(2, 2) = -0.6;
    quadratic_form(0, 1) = 0.2;
    quadratic_form(1, 0) = 0.2;
    return PolynomialField::quadratic(1.0, Vector3(0.3, -0.2, 0.1), quadratic_form);
}

std::vector<double> test_multipliers(size_t count) {
    std::vector<double> multipliers(count);

    for (size_t i = 0; i < count; ++i) {
        multipliers[i] = 0.1 * std::sin(1.7 * static_cast<double>(i));
    }

    return multipliers;
}

Vector3 finite_difference_gradient(
    const CapacityConstrainedLagrangian<>& lagrangian,
    std::vector<VectorS2> sites,
    size_t index,
    const std::vector<double>& multipliers,
    double penalty
) {
    constexpr double STEP = 1e-6;
    VectorS2 site = sites[index];
    Vector3 u = site.unitOrthogonal();
    Vector3 v = site.cross(u);
    Vector3 gradient = Vector3::Zero();

    for (const Vector3& direction : {u, v}) {
        sites[index] = (site + STEP * direction).normalized();
        double forward = lagrangian.evaluate(*build_sphere(sites), multipliers, penalty).value;

        sites[index] = (site - STEP * direction).normalized();
        double backward = lagrangian.evaluate(*build_sphere(sites), multipliers, penalty).value;

        gradient += (forward - backward) / (2.0 * STEP) * direction;
    }

    return gradient;
}

Vector3 tangential(const Vector3& gradient, const VectorS2& site) {
    return gradient - gradient.dot(site) * site;
}

}

TEST(CapacityConstrainedLagrangianTest, CapacityErrorsSumToZeroWhenTargetIsAverage) {
    auto sites = fibonacci_sites(8);
    PolynomialField field = test_field();
    CapacityConstrainedLagrangian<> lagrangian(field, field.total_mass() / sites.size());

    auto evaluation = lagrangian.evaluate(*build_sphere(sites), std::vector<double>(8, 0.0), 1.0);

    double sum = 0.0;
    for (double error : evaluation.capacity_errors) {
        sum += error;
    }

    EXPECT_NEAR(sum, 0.0, 1e-10);
}

TEST(CapacityConstrainedLagrangianTest, ValueDecomposesIntoCvtAndConstraintTerms) {
    auto sites = fibonacci_sites(8);
    PolynomialField field = test_field();
    CapacityConstrainedLagrangian<> lagrangian(field, field.total_mass() / sites.size());
    auto multipliers = test_multipliers(8);
    double penalty = 3.0;

    auto evaluation = lagrangian.evaluate(*build_sphere(sites), multipliers, penalty);

    double expected = evaluation.cvt_energy;
    for (size_t i = 0; i < 8; ++i) {
        double error = evaluation.capacity_errors[i];
        expected += multipliers[i] * error + 0.5 * penalty * error * error;
    }

    EXPECT_NEAR(evaluation.value, expected, 1e-12);
}

TEST(CapacityConstrainedLagrangianTest, GradientMatchesFiniteDifferences) {
    auto sites = fibonacci_sites(10);
    PolynomialField field = test_field();
    CapacityConstrainedLagrangian<> lagrangian(field, field.total_mass() / sites.size());
    auto multipliers = test_multipliers(10);
    double penalty = 2.5;

    auto evaluation = lagrangian.evaluate(*build_sphere(sites), multipliers, penalty);

    for (size_t index : {size_t(0), size_t(4), size_t(9)}) {
        Vector3 analytic = tangential(evaluation.site_gradients[index], sites[index]);
        Vector3 numeric = finite_difference_gradient(lagrangian, sites, index, multipliers, penalty);

        EXPECT_NEAR((analytic - numeric).norm(), 0.0, 1e-6 * std::max(1.0, numeric.norm()));
    }
}

TEST(CapacityConstrainedLagrangianTest, CvtGradientAloneMatchesFiniteDifferences) {
    auto sites = fibonacci_sites(10);
    PolynomialField field = test_field();
    CapacityConstrainedLagrangian<> lagrangian(field, field.total_mass() / sites.size());
    std::vector<double> no_multipliers(10, 0.0);

    auto evaluation = lagrangian.evaluate(*build_sphere(sites), no_multipliers, 0.0);

    Vector3 analytic = tangential(evaluation.site_gradients[3], sites[3]);
    Vector3 numeric = finite_difference_gradient(lagrangian, sites, 3, no_multipliers, 0.0);

    EXPECT_NEAR((analytic - numeric).norm(), 0.0, 1e-6 * std::max(1.0, numeric.norm()));
}

TEST(CapacityConstrainedLagrangianTest, CvtGradientIsMinusTwiceTheWeightedFirstMoment) {
    auto sites = fibonacci_sites(10);
    PolynomialField field = test_field();
    CapacityConstrainedLagrangian<> lagrangian(field, field.total_mass() / sites.size());
    auto sphere = build_sphere(sites);

    auto evaluation = lagrangian.evaluate(*sphere, std::vector<double>(10, 0.0), 0.0);
    auto states = lagrangian.cell_states(*sphere);

    for (size_t i = 0; i < sites.size(); ++i) {
        EXPECT_NEAR((evaluation.site_gradients[i] + 2.0 * states[i].first_moment).norm(), 0.0, 1e-12);
    }
}
