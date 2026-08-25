#include "voronoi_sphere_wireframe_builder.hpp"
#include "../../voronoi/spherical/core/sphere.hpp"
#include <gtest/gtest.h>
#include <cmath>

namespace geometry_art::io::mesh {
namespace {

using voronoi::spherical::Sphere;

Sphere create_small_sphere() {
    Sphere sphere;
    // Add 6 points for a simple Voronoi diagram
    sphere.insert(cgal::Point3(1, 0, 0));
    sphere.insert(cgal::Point3(-1, 0, 0));
    sphere.insert(cgal::Point3(0, 1, 0));
    sphere.insert(cgal::Point3(0, -1, 0));
    sphere.insert(cgal::Point3(0, 0, 1));
    sphere.insert(cgal::Point3(0, 0, -1));
    return sphere;
}

TEST(VoronoiSphereWireframeBuilderTest, ProducesValidMesh) {
    Sphere sphere = create_small_sphere();

    VoronoiSphereWireframeBuilder builder(0.02, 0.1);
    SurfaceMesh mesh = builder.build(sphere);

    EXPECT_GT(mesh.number_of_vertices(), 0);
    EXPECT_GT(mesh.number_of_faces(), 0);
}

TEST(VoronoiSphereWireframeBuilderTest, SmallerEdgeLengthProducesMoreVertices) {
    Sphere sphere = create_small_sphere();

    VoronoiSphereWireframeBuilder coarse_builder(0.02, 0.5);
    SurfaceMesh coarse_mesh = coarse_builder.build(sphere);

    VoronoiSphereWireframeBuilder fine_builder(0.02, 0.05);
    SurfaceMesh fine_mesh = fine_builder.build(sphere);

    EXPECT_LT(coarse_mesh.number_of_vertices(), fine_mesh.number_of_vertices());
    EXPECT_LT(coarse_mesh.number_of_faces(), fine_mesh.number_of_faces());
}

TEST(VoronoiSphereWireframeBuilderTest, VerticesAreOnCorrectRadii) {
    double thickness = 0.02;
    double half_thickness = thickness / 2.0;

    Sphere sphere = create_small_sphere();

    VoronoiSphereWireframeBuilder builder(thickness, 0.1);
    SurfaceMesh mesh = builder.build(sphere);

    double min_radius = std::numeric_limits<double>::max();
    double max_radius = std::numeric_limits<double>::lowest();

    for (auto v : mesh.vertices()) {
        auto pt = mesh.point(v);
        double radius = std::sqrt(pt.x() * pt.x() + pt.y() * pt.y() + pt.z() * pt.z());
        min_radius = std::min(min_radius, radius);
        max_radius = std::max(max_radius, radius);
    }

    EXPECT_NEAR(min_radius, 1.0 - half_thickness, 0.01);
    EXPECT_NEAR(max_radius, 1.0 + half_thickness, 0.01);
}

TEST(VoronoiSphereWireframeBuilderTest, FaceCountScalesWithArcsAndVertices) {
    Sphere small_sphere = create_small_sphere();

    Sphere larger_sphere;
    for (int i = 0; i < 20; ++i) {
        double theta = M_PI * (0.1 + 0.8 * i / 20.0);
        double phi = 2 * M_PI * i / 20.0 * 1.618;
        larger_sphere.insert(cgal::Point3(
            std::sin(theta) * std::cos(phi),
            std::sin(theta) * std::sin(phi),
            std::cos(theta)
        ));
    }

    VoronoiSphereWireframeBuilder builder(0.02, 0.1);

    SurfaceMesh small_mesh = builder.build(small_sphere);
    SurfaceMesh larger_mesh = builder.build(larger_sphere);

    EXPECT_LT(small_mesh.number_of_faces(), larger_mesh.number_of_faces());
}

} // namespace
} // namespace geometry_art::io::mesh
