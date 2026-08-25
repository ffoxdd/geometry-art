#include "triangle_mesh.hpp"
#include "polygon/polygon.hpp"
#include "arc.hpp"
#include <gtest/gtest.h>
#include <cmath>

using namespace geometry_art;
using geometry::spherical::TriangleMesh;

TEST(TriangleMeshTest, IcosahedronHasTwelveVerticesAndTwentyTriangles) {
    TriangleMesh mesh = TriangleMesh::icosphere(0);

    EXPECT_EQ(mesh.vertices.size(), 12u);
    EXPECT_EQ(mesh.triangles.size(), 20u);
}

TEST(TriangleMeshTest, SubdivisionQuadruplesTriangles) {
    TriangleMesh mesh = TriangleMesh::icosphere(2);

    EXPECT_EQ(mesh.triangles.size(), 320u);
    EXPECT_EQ(mesh.vertices.size(), 162u);
}

TEST(TriangleMeshTest, VerticesLieOnUnitSphere) {
    TriangleMesh mesh = TriangleMesh::icosphere(1);

    for (const VectorS2& vertex : mesh.vertices) {
        EXPECT_NEAR(vertex.norm(), 1.0, 1e-12);
    }
}

TEST(TriangleMeshTest, TrianglesAreCounterclockwiseAndTileTheSphere) {
    TriangleMesh mesh = TriangleMesh::icosphere(1);
    double total_area = 0.0;

    for (const auto& [a, b, c] : mesh.triangles) {
        Polygon triangle(std::vector<Arc>{
            Arc(mesh.vertices[a], mesh.vertices[b]),
            Arc(mesh.vertices[b], mesh.vertices[c]),
            Arc(mesh.vertices[c], mesh.vertices[a])
        });
        double area = triangle.area();
        EXPECT_GT(area, 0.0);
        EXPECT_LT(area, 2.0 * M_PI);
        total_area += area;
    }

    EXPECT_NEAR(total_area, 4.0 * M_PI, 1e-10);
}
