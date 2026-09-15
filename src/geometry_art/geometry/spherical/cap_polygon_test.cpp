#include "cap_polygon.hpp"
#include "polygon/polygon.hpp"
#include "../../types.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <vector>

namespace geometry_art::geometry::spherical {
namespace {

constexpr double TOLERANCE = 1e-9;

const VectorS2 X(1.0, 0.0, 0.0);
const VectorS2 Y(0.0, 1.0, 0.0);
const VectorS2 Z(0.0, 0.0, 1.0);

// The quarter of the sphere with x and y both positive: its boundary runs
// from the south pole up through +y to the north pole, then down through
// +x, counter-clockwise about the inward axis of each hemisphere.
CapPolygon lune() {
    return CapPolygon(std::vector<CapPolygon::Edge>{
        {-Z, Cap::hemisphere(X)},
        {Z, Cap::hemisphere(Y)},
    });
}

void expect_point(const VectorS2& actual, const VectorS2& expected) {
    EXPECT_NEAR((actual - expected).norm(), 0.0, TOLERANCE) <<
        "got " << actual.transpose() << " expected " << expected.transpose();
}

TEST(CapPolygonTest, OfAGreatCirclePolygonKeepsItsVerticesOnHemispheres) {
    std::vector<VectorS2> corners{X, Y, Z};
    CapPolygon region = CapPolygon::of(*polygon::Polygon::from_points(corners));

    ASSERT_EQ(region.size(), 3u);

    for (const CapPolygon::Edge& edge : region.edges()) {
        EXPECT_DOUBLE_EQ(edge.cap.offset, 0.0);
        EXPECT_NEAR(edge.cap.axis.norm(), 1.0, TOLERANCE);
    }
}

TEST(CapPolygonTest, ACapContainingEverythingChangesNothing) {
    CapPolygon region = lune();
    region.clip(Cap::hemisphere(X + Y));

    ASSERT_EQ(region.size(), 2u);
    expect_point(region.edges()[0].source, -Z);
    expect_point(region.edges()[1].source, Z);
}

TEST(CapPolygonTest, ACapMissingEverythingEmptiesIt) {
    CapPolygon region = lune();
    region.clip(Cap::hemisphere_inset_by(-X - Y, 0.1));

    EXPECT_TRUE(region.empty());
}

TEST(CapPolygonTest, AnInsetCapMovesAnEdgeOntoTheSmallCircle) {
    double inset = 0.3;
    CapPolygon region = lune();
    region.clip(Cap::hemisphere_inset_by(X, inset));

    ASSERT_EQ(region.size(), 2u);
    expect_point(region.edges()[0].source, VectorS2(std::sin(inset), 0.0, std::cos(inset)));
    expect_point(region.edges()[1].source, VectorS2(std::sin(inset), 0.0, -std::cos(inset)));
    EXPECT_DOUBLE_EQ(region.edges()[1].cap.offset, std::sin(inset));
    expect_point(region.point(1, 0.5), VectorS2(std::sin(inset), std::cos(inset), 0.0));
    EXPECT_NEAR(region.length(1), M_PI * std::cos(inset), TOLERANCE);
}

TEST(CapPolygonTest, AnEdgeNickedInItsMiddleKeepsBothEnds) {
    double reach = 20.0 * M_PI / 180.0;
    CapPolygon region = lune();
    region.clip(Cap{-X, -std::cos(reach)});

    ASSERT_EQ(region.size(), 4u);
    expect_point(region.edges()[0].source, -Z);
    expect_point(region.edges()[1].source, Z);
    expect_point(region.edges()[2].source, VectorS2(std::cos(reach), 0.0, std::sin(reach)));
    expect_point(region.edges()[3].source, VectorS2(std::cos(reach), 0.0, -std::sin(reach)));
    expect_point(region.point(2, 0.5), VectorS2(std::cos(reach), std::sin(reach), 0.0));
}

TEST(CapPolygonTest, APointOnAnEdgeAdvancesCounterClockwiseAboutItsAxis) {
    CapPolygon region = lune();

    expect_point(region.point(0, 0.5), Y);
    expect_point(region.point(1, 0.5), X);
    EXPECT_NEAR(region.length(0), M_PI, TOLERANCE);
}

} // namespace
} // namespace geometry_art::geometry::spherical
