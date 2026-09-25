#include "welder.hpp"
#include "../types.hpp"
#include <gtest/gtest.h>

namespace geometry_art::skeleton {
namespace {

TEST(WelderTest, ThePositionTwiceIsOneVertex) {
    Welder welder(1.0);

    EXPECT_EQ(welder.vertex(Vector3(0.1, 0.2, 0.3)), welder.vertex(Vector3(0.1, 0.2, 0.3)));
}

TEST(WelderTest, PositionsEitherSideOfARoundingBoundaryWeld) {
    Welder welder(1.0);

    EXPECT_EQ(welder.vertex(Vector3(0.7999999995, 0.6, 0.0)), welder.vertex(Vector3(0.79999999949999989, 0.6, 0.0)));
    EXPECT_EQ(welder.vertex(Vector3(0.0000004, 0.0, 0.0)), welder.vertex(Vector3(0.0000011, 0.0, 0.0)));
}

TEST(WelderTest, PositionsTooCloseForTheOutputWeld) {
    Welder welder(1.0);

    EXPECT_EQ(welder.vertex(Vector3(0.6, 0.8, 0.0)), welder.vertex(Vector3(0.600000003, 0.800000003, 0.0)));
}

TEST(WelderTest, DistinctPositionsStayApart) {
    Welder welder(1.0);

    EXPECT_NE(welder.vertex(Vector3(0.6, 0.8, 0.0)), welder.vertex(Vector3(0.6001, 0.8, 0.0)));
}

TEST(WelderTest, ATriangleWithAWeldedEdgeIsDropped) {
    Welder welder(1.0);
    VertexIndex a = welder.vertex(Vector3(0.0, 0.0, 0.0));
    VertexIndex b = welder.vertex(Vector3(0.0000001, 0.0, 0.0));
    VertexIndex c = welder.vertex(Vector3(0.0, 1.0, 0.0));

    welder.triangle(a, b, c);

    EXPECT_EQ(welder.take().number_of_faces(), 0u);
}

TEST(WelderTest, PositionsAreStoredAtTheOutputScale) {
    Welder welder(10.0);
    VertexIndex index = welder.vertex(Vector3(0.1, 0.2, 0.3));
    SurfaceMesh mesh = welder.take();

    EXPECT_NEAR(mesh.point(index).x(), 1.0, 1e-12);
    EXPECT_NEAR(mesh.point(index).z(), 3.0, 1e-12);
}

} // namespace
} // namespace geometry_art::skeleton
