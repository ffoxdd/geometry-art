#include "annulus.hpp"
#include "../types.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <vector>

namespace geometry_art::skeleton {
namespace {

std::vector<Vector2> square(double half, const Vector2& center = Vector2::Zero()) {
    return {
        center + Vector2(-half, -half),
        center + Vector2(half, -half),
        center + Vector2(half, half),
        center + Vector2(-half, half),
    };
}

double signed_area(const Vector2& a, const Vector2& b, const Vector2& c) {
    return 0.5 * ((b.x() - a.x()) * (c.y() - a.y()) - (b.y() - a.y()) * (c.x() - a.x()));
}

double total_area(const std::vector<Vector2>& outer, const std::vector<Vector2>& inner, const std::vector<Triangle>& triangles) {
    std::vector<Vector2> points = outer;
    points.insert(points.end(), inner.begin(), inner.end());
    double total = 0.0;

    for (const Triangle& triangle : triangles) {
        double area = signed_area(points[triangle[0]], points[triangle[1]], points[triangle[2]]);
        EXPECT_GT(area, 0.0);
        total += area;
    }

    return total;
}

TEST(AnnulusTest, NestedSquaresAreCoveredByEightCounterClockwiseTriangles) {
    std::vector<Vector2> outer = square(2.0);
    std::vector<Vector2> inner = square(1.0);

    std::vector<Triangle> triangles = triangulate_annulus(outer, inner);

    EXPECT_EQ(triangles.size(), 8u);
    EXPECT_NEAR(total_area(outer, inner, triangles), 16.0 - 4.0, 1e-12);
}

TEST(AnnulusTest, AnOffCentreInnerLoopStillCoversTheRing) {
    std::vector<Vector2> outer = square(2.0);
    std::vector<Vector2> inner = square(0.5, Vector2(1.0, -0.8));

    std::vector<Triangle> triangles = triangulate_annulus(outer, inner);

    EXPECT_EQ(triangles.size(), 8u);
    EXPECT_NEAR(total_area(outer, inner, triangles), 16.0 - 1.0, 1e-12);
}

TEST(AnnulusTest, LoopsOfDifferentLengthsZipTogether) {
    std::vector<Vector2> outer;

    for (int step = 0; step < 12; ++step) {
        double angle = TWO_PI * step / 12.0;
        outer.emplace_back(3.0 * std::cos(angle), 3.0 * std::sin(angle));
    }

    std::vector<Vector2> inner = square(1.0);
    std::vector<Triangle> triangles = triangulate_annulus(outer, inner);

    EXPECT_EQ(triangles.size(), 16u);
    EXPECT_NEAR(total_area(outer, inner, triangles), 0.5 * 12.0 * 9.0 * std::sin(TWO_PI / 12.0) - 4.0, 1e-12);
}

TEST(AnnulusTest, WithoutAnInnerLoopTheOuterOneIsFanned) {
    std::vector<Vector2> outer = square(1.0);

    std::vector<Triangle> triangles = triangulate_annulus(outer, {});

    EXPECT_EQ(triangles.size(), 2u);
    EXPECT_NEAR(total_area(outer, {}, triangles), 4.0, 1e-12);
}

TEST(AnnulusTest, TooFewOuterVerticesGiveNothing) {
    EXPECT_TRUE(triangulate_annulus({Vector2(0.0, 0.0), Vector2(1.0, 0.0)}, {}).empty());
}

} // namespace
} // namespace geometry_art::skeleton
