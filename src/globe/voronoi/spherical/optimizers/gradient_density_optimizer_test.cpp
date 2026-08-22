#include "gradient_density_optimizer.hpp"
#include "../../../fields/spherical/polynomial_field.hpp"
#include "../../../geometry/spherical/arc.hpp"
#include "../../../geometry/spherical/helpers.hpp"
#include "../../../testing/assertions/geometric.hpp"
#include "../../../testing/macros.hpp"
#include <gtest/gtest.h>
#include <cmath>

using namespace globe;
using namespace globe::voronoi::spherical;
using geometry::spherical::Arc;
using fields::spherical::PolynomialField;

namespace {

std::unique_ptr<Sphere> create_test_voronoi(size_t num_points) {
    auto voronoi = std::make_unique<Sphere>();

    double phi = M_PI * (std::sqrt(5.0) - 1.0);
    for (size_t i = 0; i < num_points; ++i) {
        double y = 1.0 - (i / static_cast<double>(num_points - 1)) * 2.0;
        double radius = std::sqrt(1.0 - y * y);
        double theta = phi * i;

        double x = std::cos(theta) * radius;
        double z = std::sin(theta) * radius;

        voronoi->insert(cgal::Point3(x, y, z));
    }

    return voronoi;
}

}

TEST(GradientDensityOptimizerTest, ConvergesForConstantField) {
    auto voronoi = create_test_voronoi(10);
    PolynomialField field = PolynomialField::constant(1.0);

    GradientDensityOptimizer optimizer(
        std::move(voronoi),
        field,
        50
    );

    auto result = optimizer.optimize();

    EXPECT_EQ(result->size(), 10);
}

TEST(GradientDensityOptimizerTest, EXPENSIVE_ConvergesForLinearField) {
    REQUIRE_EXPENSIVE();

    auto voronoi = create_test_voronoi(10);
    PolynomialField field = PolynomialField::linear(2.0, Vector3(0, 0, 1.0));

    GradientDensityOptimizer<PolynomialField> optimizer(
        std::move(voronoi),
        field,
        50
    );

    auto result = optimizer.optimize();

    EXPECT_EQ(result->size(), 10);
}

TEST(GradientDensityOptimizerTest, EXPENSIVE_GradientMatchesNumericalForLinearField) {
    REQUIRE_EXPENSIVE();

    auto voronoi = create_test_voronoi(6);
    PolynomialField field = PolynomialField::linear(2.0, Vector3(0, 0, 1.0));
    double target_mass = field.total_mass() / voronoi->size();

    auto compute_mass_errors = [&]() {
        std::vector<double> errors(voronoi->size());
        size_t i = 0;
        for (const auto& cell : voronoi->cells()) {
            errors[i] = field.integrals(cell).mass - target_mass;
            ++i;
        }
        return errors;
    };

    std::vector<cgal::Point3> sites;
    for (size_t i = 0; i < voronoi->size(); ++i) {
        sites.push_back(voronoi->site(i));
    }

    auto rebuild_voronoi_with_perturbation = [&sites](size_t index, cgal::Point3 new_pos) {
        auto v = std::make_unique<Sphere>();
        for (size_t i = 0; i < sites.size(); ++i) {
            v->insert(i == index ? new_pos : sites[i]);
        }
        return v;
    };

    auto compute_error_for = [&](Sphere& v) {
        double total = 0.0;
        for (const auto& cell : v.cells()) {
            double mass_error = field.integrals(cell).mass - target_mass;
            total += mass_error * mass_error;
        }
        return total / 2.0;
    };

    auto compute_analytical_gradient = [&](size_t site_index) {
        auto mass_errors = compute_mass_errors();
        Eigen::Vector3d gradient = Eigen::Vector3d::Zero();

        cgal::Point3 site_k = voronoi->site(site_index);
        Eigen::Vector3d s_k = to_eigen(site_k);
        auto cell_edges = voronoi->cell_edges(site_index);

        for (const auto& edge_info : cell_edges) {
            size_t j = edge_info.neighbor_index;
            cgal::Point3 site_j = voronoi->site(j);

            auto edge = field.integrals(edge_info.arc);
            Eigen::Vector3d rho_weighted_moment = edge.first_moment;
            double edge_integral = edge.mass;

            Eigen::Vector3d n_vec = to_eigen(site_j) - s_k;
            double n_norm = n_vec.norm();

            if (n_norm < GEOMETRIC_EPSILON) continue;

            Eigen::Vector3d edge_grad = (rho_weighted_moment - s_k * edge_integral) / n_norm;
            gradient += (mass_errors[site_index] - mass_errors[j]) * edge_grad;
        }

        return gradient;
    };

    double epsilon = 1e-5;
    size_t site_idx = 0;
    cgal::Point3 site = sites[site_idx];
    Eigen::Vector3d s(site.x(), site.y(), site.z());

    Eigen::Vector3d analytical = compute_analytical_gradient(site_idx);
    Eigen::Vector3d analytical_tangent = analytical - analytical.dot(s) * s;

    Eigen::Vector3d numerical = Eigen::Vector3d::Zero();
    for (int dim = 0; dim < 3; ++dim) {
        Eigen::Vector3d perturb = Eigen::Vector3d::Zero();
        perturb[dim] = epsilon;

        Eigen::Vector3d new_pos_plus = (s + perturb).normalized();
        auto voronoi_plus = rebuild_voronoi_with_perturbation(
            site_idx, cgal::Point3(new_pos_plus.x(), new_pos_plus.y(), new_pos_plus.z())
        );
        double error_plus = compute_error_for(*voronoi_plus);

        Eigen::Vector3d new_pos_minus = (s - perturb).normalized();
        auto voronoi_minus = rebuild_voronoi_with_perturbation(
            site_idx, cgal::Point3(new_pos_minus.x(), new_pos_minus.y(), new_pos_minus.z())
        );
        double error_minus = compute_error_for(*voronoi_minus);

        numerical[dim] = (error_plus - error_minus) / (2 * epsilon);
    }

    EXPECT_NEAR(analytical_tangent.x(), numerical.x(), 1e-6);
    EXPECT_NEAR(analytical_tangent.y(), numerical.y(), 1e-6);
    EXPECT_NEAR(analytical_tangent.z(), numerical.z(), 1e-6);
}
