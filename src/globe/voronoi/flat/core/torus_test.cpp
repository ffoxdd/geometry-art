#include "torus.hpp"
#include "../../../testing/macros.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <cstddef>
#include <set>
#include <vector>

using namespace globe;
using namespace globe::voronoi::flat;

namespace {

// Low-discrepancy sites covering the rectangle, so no cell is degenerate.
std::vector<Vector2> scattered_sites(size_t count, double width, double height) {
    std::vector<Vector2> sites;
    double golden = 0.6180339887498949;

    for (size_t k = 0; k < count; ++k) {
        double x = std::fmod(0.13 + golden * static_cast<double>(k), 1.0) * width;
        double y = (static_cast<double>(k) + 0.5) / static_cast<double>(count) * height;
        sites.emplace_back(x, y);
    }

    return sites;
}

Torus scattered_torus(size_t count, double width, double height) {
    Torus torus(width, height);

    for (const Vector2& site : scattered_sites(count, width, height)) {
        torus.insert(site);
    }

    return torus;
}

} // namespace

TEST(TorusTest, CanonicalWrapsIntoTheRectangle) {
    Torus torus(2.0, 1.0);

    EXPECT_NEAR(torus.canonical(Vector2(2.5, -0.25)).x(), 0.5, 1e-12);
    EXPECT_NEAR(torus.canonical(Vector2(2.5, -0.25)).y(), 0.75, 1e-12);
    EXPECT_NEAR(torus.canonical(Vector2(-0.5, 1.0)).x(), 1.5, 1e-12);
    EXPECT_NEAR(torus.canonical(Vector2(-0.5, 1.0)).y(), 0.0, 1e-12);
}

// A regular grid of sites tessellates the torus into identical rectangles,
// the one case where every cell is known in closed form.
TEST(TorusTest, GridSitesGiveRectangularCellsOfEqualArea) {
    Torus torus(2.0, 1.0);

    for (int row = 0; row < 2; ++row) {
        for (int column = 0; column < 4; ++column) {
            torus.insert(Vector2(0.25 + 0.5 * column, 0.25 + 0.5 * row));
        }
    }

    for (size_t index = 0; index < torus.size(); ++index) {
        EXPECT_NEAR(torus.cell(index).area(), 0.25, 1e-9);
    }
}

TEST(TorusTest, CellAreasPartitionTheRectangle) {
    Torus torus = scattered_torus(17, 2.0, 1.0);
    double total = 0.0;

    for (size_t index = 0; index < torus.size(); ++index) {
        total += torus.cell(index).area();
    }

    EXPECT_NEAR(total, torus.area(), 1e-9);
}

TEST(TorusTest, NeighborsAreSymmetric) {
    Torus torus = scattered_torus(13, 1.5, 1.0);

    for (size_t index = 0; index < torus.size(); ++index) {
        for (const CellEdgeInfo& edge : torus.cell_edges(index)) {
            bool found = false;

            for (const CellEdgeInfo& back : torus.cell_edges(edge.neighbor_index)) {
                found = found || back.neighbor_index == index;
            }

            EXPECT_TRUE(found) << index << " -> " << edge.neighbor_index;
        }
    }
}

// The neighbor's position is reported in the cell's own chart: exactly one
// period image of the neighbor's canonical site, and the closer one.
TEST(TorusTest, NeighborPositionsAreChartConsistent) {
    Torus torus = scattered_torus(13, 1.5, 1.0);

    for (size_t index = 0; index < torus.size(); ++index) {
        for (const CellEdgeInfo& edge : torus.cell_edges(index)) {
            Vector2 canonical = torus.site(edge.neighbor_index);
            Vector2 offset = edge.neighbor_position - canonical;
            double x_periods = offset.x() / torus.width();
            double y_periods = offset.y() / torus.height();

            EXPECT_NEAR(x_periods, std::round(x_periods), 1e-9);
            EXPECT_NEAR(y_periods, std::round(y_periods), 1e-9);
        }
    }
}

// A bisector's endpoints are equidistant from the two sites, in the chart
// where all three were reported.
TEST(TorusTest, BoundariesAreEquidistantFromBothSites) {
    Torus torus = scattered_torus(11, 1.0, 1.0);

    for (size_t index = 0; index < torus.size(); ++index) {
        Vector2 own = torus.site(index);

        for (const CellEdgeInfo& edge : torus.cell_edges(index)) {
            for (const Vector2& endpoint : {edge.boundary.source(), edge.boundary.target()}) {
                double to_own = (endpoint - own).norm();
                double to_neighbor = (endpoint - edge.neighbor_position).norm();

                EXPECT_NEAR(to_own, to_neighbor, 1e-9);
            }
        }
    }
}

// With this few sites a cell wraps far enough to border itself: the seam
// where it meets its own period image is a genuine equidistance locus, and
// it must appear in the edge list -- twice, once per direction, so the two
// cuts' sweep contributions cancel exactly.
TEST(TorusTest, ACellCanBorderItselfAcrossTheSeam) {
    Torus torus(1.0, 1.0);
    torus.insert(Vector2(0.25, 0.5));
    torus.insert(Vector2(0.75, 0.5));

    EXPECT_NEAR(torus.cell(0).area(), 0.5, 1e-9);
    EXPECT_NEAR(torus.cell(1).area(), 0.5, 1e-9);

    std::set<size_t> neighbors;
    size_t self_edges = 0;

    for (const CellEdgeInfo& edge : torus.cell_edges(0)) {
        neighbors.insert(edge.neighbor_index);
        self_edges += edge.neighbor_index == 0 ? 1 : 0;
    }

    EXPECT_EQ(neighbors, (std::set<size_t>{0, 1}));
    EXPECT_EQ(self_edges % 2, 0u);
}

TEST(TorusTest, RebuiltPreservesSitesUpToWrapping) {
    Torus torus = scattered_torus(9, 2.0, 1.0);
    std::vector<Vector3> moved;

    for (size_t index = 0; index < torus.size(); ++index) {
        moved.push_back(torus.site_vector(index) + Vector3(2.0, -1.0, 0.0));
    }

    auto rebuilt = torus.rebuilt(moved);

    ASSERT_EQ(rebuilt->size(), torus.size());

    for (size_t index = 0; index < torus.size(); ++index) {
        EXPECT_NEAR((rebuilt->site(index) - torus.site(index)).norm(), 0.0, 1e-12);
    }
}
