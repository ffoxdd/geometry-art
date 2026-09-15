#include "builder.hpp"
#include "flat_outliner.hpp"
#include "spherical_outliner.hpp"
#include "../cgal/types.hpp"
#include "../geometry/planar/domain.hpp"
#include "../geometry/planar/polygon.hpp"
#include "../io/mesh/types.hpp"
#include "../testing/flat_scatter.hpp"
#include "../testing/macros.hpp"
#include "../types.hpp"
#include "../voronoi/flat/core/diagram.hpp"
#include "../voronoi/spherical/core/sphere.hpp"
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <gtest/gtest.h>
#include <cmath>
#include <memory>
#include <vector>

namespace geometry_art::skeleton {
namespace {

namespace PolygonMeshProcessing = CGAL::Polygon_mesh_processing;

using geometry::planar::Domain;
using geometry_art::testing::flat_scatter;
using io::mesh::SurfaceMesh;
using voronoi::flat::CellClipper;
using voronoi::flat::Diagram;
using voronoi::spherical::Sphere;

const Parameters FINE{0.1, 0.05, 0.2, 1.0};

std::vector<Vector2> grid_sites() {
    std::vector<Vector2> sites;

    for (int row = 0; row < 2; ++row) {
        for (int column = 0; column < 4; ++column) {
            sites.emplace_back(0.25 + 0.5 * column, 0.25 + 0.5 * row);
        }
    }

    return sites;
}

Sphere fibonacci_sphere(int count) {
    Sphere sphere;
    double golden = M_PI * (3.0 - std::sqrt(5.0));

    for (int index = 0; index < count; ++index) {
        double z = 1.0 - 2.0 * (index + 0.5) / count;
        double ring = std::sqrt(1.0 - z * z);
        double angle = golden * index;
        sphere.insert(cgal::Point3(ring * std::cos(angle), ring * std::sin(angle), z));
    }

    return sphere;
}

double region_area(const Diagram& diagram, double half_width) {
    const Domain& domain = diagram.domain();
    double area = (domain.width + (domain.wrapped(0) ? 0.0 : 2.0 * half_width)) *
        (domain.height + (domain.wrapped(1) ? 0.0 : 2.0 * half_width));

    for (size_t index = 0; index < diagram.size(); ++index) {
        CellClipper inset = diagram.offset_cell(index, half_width, half_width);
        std::vector<Vector2> vertices;

        for (const CellClipper::Edge& edge : inset.edges()) {
            vertices.push_back(edge.source);
        }

        if (vertices.size() >= 3) {
            area -= geometry::planar::Polygon(vertices).area();
        }
    }

    return area;
}

void expect_closed(const SurfaceMesh& mesh) {
    EXPECT_GT(mesh.number_of_faces(), 0u);
    EXPECT_TRUE(CGAL::is_closed(mesh));
    EXPECT_TRUE(CGAL::is_valid_polygon_mesh(mesh));
}

TEST(BuilderTest, ThePlaneSkeletonIsAClosedSolidOfTheRegionsVolume) {
    Diagram diagram(Domain::plane(2.0, 1.0), grid_sites());
    SurfaceMesh mesh = Builder(FINE).build(FlatOutliner(diagram));

    expect_closed(mesh);
    EXPECT_NEAR(PolygonMeshProcessing::volume(mesh), FINE.bar_thickness * region_area(diagram, 0.5 * FINE.bar_width), 1e-9);
}

TEST(BuilderTest, TheScaleMultipliesEveryPosition) {
    Diagram diagram(Domain::plane(2.0, 1.0), grid_sites());
    Parameters scaled = FINE;
    scaled.scale = 10.0;
    SurfaceMesh mesh = Builder(scaled).build(FlatOutliner(diagram));

    expect_closed(mesh);
    EXPECT_NEAR(PolygonMeshProcessing::volume(mesh), 1000.0 * FINE.bar_thickness * region_area(diagram, 0.5 * FINE.bar_width), 1e-6);
}

TEST(BuilderTest, ACellTheBarsFillBecomesASolidBlock) {
    Diagram diagram(Domain::plane(0.2, 0.2), {Vector2(0.1, 0.1)});
    Parameters wide{0.5, 0.1, 0.0, 1.0};
    SurfaceMesh mesh = Builder(wide).build(FlatOutliner(diagram));

    expect_closed(mesh);
    EXPECT_NEAR(PolygonMeshProcessing::volume(mesh), 0.1 * 0.7 * 0.7, 1e-12);
}

TEST(BuilderTest, TheCylinderSkeletonClosesAcrossTheSeamAndAlongTheRims) {
    Diagram diagram(Domain::cylinder(2.0, 1.0), grid_sites());
    SurfaceMesh mesh = Builder(FINE).build(FlatOutliner(diagram));

    expect_closed(mesh);
    EXPECT_GT(PolygonMeshProcessing::volume(mesh), 0.0);
}

TEST(BuilderTest, TheTorusSkeletonClosesAcrossBothSeams) {
    Diagram diagram(Domain::torus(2.0, 1.0), grid_sites());
    SurfaceMesh mesh = Builder(FINE).build(FlatOutliner(diagram));

    expect_closed(mesh);
    EXPECT_GT(PolygonMeshProcessing::volume(mesh), 0.0);
}

TEST(BuilderTest, TheSphereSkeletonLiesBetweenItsTwoRadii) {
    Sphere sphere = fibonacci_sphere(12);
    SurfaceMesh mesh = Builder(FINE).build(SphericalOutliner(sphere));

    expect_closed(mesh);
    EXPECT_GT(PolygonMeshProcessing::volume(mesh), 0.0);

    double inner = 1.0 - 0.5 * FINE.bar_thickness;
    double outer = 1.0 + 0.5 * FINE.bar_thickness;

    for (auto vertex : mesh.vertices()) {
        double radius = to_vector3(mesh.point(vertex)).norm();
        EXPECT_TRUE(std::abs(radius - inner) < 1e-9 || std::abs(radius - outer) < 1e-9) << radius;
    }
}

TEST(BuilderTest, TooFewSitesGiveAnEmptyMesh) {
    Sphere sphere;
    sphere.insert(cgal::Point3(1, 0, 0));
    sphere.insert(cgal::Point3(-1, 0, 0));

    EXPECT_EQ(Builder(FINE).build(SphericalOutliner(sphere)).number_of_faces(), 0u);
}

TEST(BuilderTest, EXPENSIVE_SkeletonsDoNotSelfIntersect) {
    REQUIRE_EXPENSIVE();

    Parameters bars{0.06, 0.03, 0.05, 1.0};

    Diagram plane(Domain::plane(2.0, 1.0), flat_scatter(40, 2.0, 1.0));
    Diagram cylinder(Domain::cylinder(2.0, 1.0), flat_scatter(40, 2.0, 1.0));
    Diagram torus(Domain::torus(2.0, 1.0), flat_scatter(40, 2.0, 1.0));
    Sphere sphere = fibonacci_sphere(40);

    for (SurfaceMesh mesh : {
        Builder(bars).build(FlatOutliner(plane)),
        Builder(bars).build(FlatOutliner(cylinder)),
        Builder(bars).build(FlatOutliner(torus)),
        Builder(bars).build(SphericalOutliner(sphere)),
    }) {
        expect_closed(mesh);
        EXPECT_FALSE(PolygonMeshProcessing::does_self_intersect(mesh));
    }
}

} // namespace
} // namespace geometry_art::skeleton
