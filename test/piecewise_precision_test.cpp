#include "geometry_art/fields/spherical/piecewise_polynomial_field.hpp"
#include "geometry_art/fields/spherical/polynomial_field.hpp"
#include "geometry_art/generators/cartesian/random_point_generator.hpp"
#include "geometry_art/generators/spherical/random_point_generator.hpp"
#include "geometry_art/geometry/cartesian/bounding_box_sampler/uniform_bounding_box_sampler.hpp"
#include "geometry_art/geometry/spherical/triangle_mesh.hpp"
#include "geometry_art/math/interval_sampler/uniform_interval_sampler.hpp"
#include "geometry_art/testing/macros.hpp"
#include "geometry_art/voronoi/spherical/core/random_builder.hpp"
#include "geometry_art/voronoi/spherical/core/sphere.hpp"
#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <memory>

using namespace geometry_art;
using fields::spherical::PiecewisePolynomialField;
using fields::spherical::PolynomialField;
using fields::RegionIntegrals;
using geometry::spherical::TriangleMesh;
using voronoi::spherical::RandomBuilder;
using voronoi::spherical::Sphere;

namespace {

using SeededBoundingBoxSampler = UniformBoundingBoxSampler<UniformIntervalSampler>;
using SeededCartesianGenerator = generators::cartesian::RandomPointGenerator<SeededBoundingBoxSampler>;
using SeededPointGenerator = generators::spherical::RandomPointGenerator<SeededCartesianGenerator>;

std::unique_ptr<Sphere> random_sphere(unsigned int seed, int count) {
    SeededPointGenerator generator{
        SeededCartesianGenerator{SeededBoundingBoxSampler{UniformIntervalSampler{seed}}},
        UniformSphericalBoundingBoxSampler<>{}
    };

    return RandomBuilder<SeededPointGenerator>(generator).build(count);
}

PolynomialField quadratic() {
    Eigen::Matrix3d matrix = Eigen::Matrix3d::Zero();
    matrix(2, 2) = -0.9;
    return PolynomialField::quadratic(1.0, Vector3::Zero(), matrix);
}

struct Deviation {
    double mass = 0.0;
    double first_moment = 0.0;
    double arc_mass = 0.0;
    double arc_first_moment = 0.0;
    double second_moment = 0.0;
};

Deviation measure(const Sphere& sphere, const PolynomialField& global, const PiecewisePolynomialField& piecewise) {
    Deviation worst;

    for (size_t k = 0; k < sphere.size(); ++k) {
        Polygon cell = sphere.cell(k);
        RegionIntegrals expected = global.integrals(cell);
        RegionIntegrals actual = piecewise.integrals(cell);
        worst.mass = std::max(worst.mass, std::abs(actual.mass - expected.mass) / expected.mass);
        worst.first_moment = std::max(worst.first_moment, (actual.first_moment - expected.first_moment).norm() / expected.mass);

        for (const auto& edge : sphere.cell_edges(k)) {
            RegionIntegrals expected_arc = global.integrals(edge.arc);
            RegionIntegrals actual_arc = piecewise.integrals(edge.arc);
            double scale = std::max(expected_arc.mass, 1e-300);
            worst.arc_mass = std::max(worst.arc_mass, std::abs(actual_arc.mass - expected_arc.mass) / scale);
            worst.arc_first_moment = std::max(worst.arc_first_moment, (actual_arc.first_moment - expected_arc.first_moment).norm() / scale);
            worst.second_moment = std::max(worst.second_moment, (piecewise.second_moment(edge.arc) - global.second_moment(edge.arc)).norm() / scale);
        }
    }

    return worst;
}

} // namespace

TEST(PiecewisePrecisionTest, EXPENSIVE_PiecewiseQuadraticMatchesGlobalOnEveryCellAndEdge) {
    REQUIRE_EXPENSIVE();
    PolynomialField global = quadratic();
    PiecewisePolynomialField piecewise = PiecewisePolynomialField::sample(TriangleMesh::icosphere(4), 2, global);

    for (unsigned int seed : {7u, 11u, 13u}) {
        auto sphere = random_sphere(seed, 200);
        Deviation worst = measure(*sphere, global, piecewise);
        EXPECT_LT(worst.mass, 1e-11) << "seed " << seed;
        EXPECT_LT(worst.first_moment, 1e-11) << "seed " << seed;
        EXPECT_LT(worst.arc_mass, 1e-11) << "seed " << seed;
        EXPECT_LT(worst.arc_first_moment, 1e-11) << "seed " << seed;
        EXPECT_LT(worst.second_moment, 1e-11) << "seed " << seed;
    }

    EXPECT_NEAR(piecewise.total_mass(), global.total_mass(), 1e-12);
}
