#include "diagram.hpp"
#include "../../../geometry/planar/domain.hpp"
#include "../../../testing/flat_scatter.hpp"
#include "../../../testing/macros.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <cstddef>
#include <set>
#include <variant>
#include <vector>

using namespace geometry_art;
using namespace geometry_art::voronoi::flat;
using geometry::planar::Domain;
using geometry_art::testing::diagram_of;
using geometry_art::testing::flat_scatter;
using geometry_art::testing::scattered_plane;
using geometry_art::testing::scattered_torus;

namespace {

std::vector<Vector2> grid_sites() {
    std::vector<Vector2> sites;

    for (int row = 0; row < 2; ++row) {
        for (int column = 0; column < 4; ++column) {
            sites.emplace_back(0.25 + 0.5 * column, 0.25 + 0.5 * row);
        }
    }

    return sites;
}

double total_cell_area(const Diagram& diagram) {
    double total = 0.0;

    for (size_t index = 0; index < diagram.size(); ++index) {
        total += diagram.cell(index).area();
    }

    return total;
}

} // namespace

TEST(FlatDiagramTest, CanonicalWrapsIntoTheRectangle) {
    Diagram diagram(Domain::torus(2.0, 1.0), {});

    EXPECT_NEAR(diagram.canonical(Vector2(2.5, -0.25)).x(), 0.5, 1e-12);
    EXPECT_NEAR(diagram.canonical(Vector2(2.5, -0.25)).y(), 0.75, 1e-12);
    EXPECT_NEAR(diagram.canonical(Vector2(-0.5, 1.0)).x(), 1.5, 1e-12);
    EXPECT_NEAR(diagram.canonical(Vector2(-0.5, 1.0)).y(), 0.0, 1e-12);
}

// A regular grid of sites tessellates the torus into identical rectangles,
// the one case where every cell is known in closed form.
TEST(FlatDiagramTest, GridSitesGiveRectangularCellsOfEqualArea) {
    Diagram diagram(Domain::torus(2.0, 1.0), grid_sites());

    for (size_t index = 0; index < diagram.size(); ++index) {
        EXPECT_NEAR(diagram.cell(index).area(), 0.25, 1e-9);
    }
}

TEST(FlatDiagramTest, CellAreasPartitionTheRectangle) {
    auto diagram = scattered_torus(17, 2.0, 1.0);

    EXPECT_NEAR(total_cell_area(*diagram), diagram->area(), 1e-9);
}

TEST(FlatDiagramTest, NeighborsAreSymmetric) {
    auto diagram = scattered_torus(13, 1.5, 1.0);

    for (size_t index = 0; index < diagram->size(); ++index) {
        for (const CellEdgeInfo& edge : diagram->cell_edges(index)) {
            bool found = false;

            for (const CellEdgeInfo& back : diagram->cell_edges(edge.neighbor_index)) {
                found = found || back.neighbor_index == index;
            }

            EXPECT_TRUE(found) << index << " -> " << edge.neighbor_index;
        }
    }
}

// The neighbor's position is reported in the cell's own chart: exactly one
// period image of the neighbor's canonical site, and the closer one.
TEST(FlatDiagramTest, NeighborPositionsAreChartConsistent) {
    auto diagram = scattered_torus(13, 1.5, 1.0);

    for (size_t index = 0; index < diagram->size(); ++index) {
        for (const CellEdgeInfo& edge : diagram->cell_edges(index)) {
            Vector2 canonical = diagram->site(edge.neighbor_index);
            Vector2 offset = edge.neighbor_position - canonical;
            double x_periods = offset.x() / diagram->width();
            double y_periods = offset.y() / diagram->height();

            EXPECT_NEAR(x_periods, std::round(x_periods), 1e-9);
            EXPECT_NEAR(y_periods, std::round(y_periods), 1e-9);
        }
    }
}

// A bisector's endpoints are equidistant from the two sites, in the chart
// where all three were reported.
TEST(FlatDiagramTest, BoundariesAreEquidistantFromBothSites) {
    auto diagram = scattered_torus(11, 1.0, 1.0);

    for (size_t index = 0; index < diagram->size(); ++index) {
        Vector2 own = diagram->site(index);

        for (const CellEdgeInfo& edge : diagram->cell_edges(index)) {
            for (const Vector2& endpoint : {edge.boundary.source(), edge.boundary.target()}) {
                double to_own = (endpoint - own).norm();
                double to_neighbor = (endpoint - edge.neighbor_position).norm();

                EXPECT_NEAR(to_own, to_neighbor, 1e-9);
            }
        }
    }
}

// The cut on the far side of an endpoint is the third constraint on it:
// where it names a site, the endpoint is equidistant from that site too.
TEST(FlatDiagramTest, EndpointCutsNameTheThirdSite) {
    auto diagram = scattered_torus(11, 1.0, 1.0);

    for (size_t index = 0; index < diagram->size(); ++index) {
        Vector2 own = diagram->site(index);

        for (const CellEdgeInfo& edge : diagram->cell_edges(index)) {
            const Bisector& source = std::get<Bisector>(edge.source_cut);
            const Bisector& target = std::get<Bisector>(edge.target_cut);

            EXPECT_NEAR((edge.boundary.source() - own).norm(), (edge.boundary.source() - source.neighbor_position).norm(), 1e-9);
            EXPECT_NEAR((edge.boundary.target() - own).norm(), (edge.boundary.target() - target.neighbor_position).norm(), 1e-9);
        }
    }
}

// With this few sites a cell wraps far enough to border itself: the seam
// where it meets its own period image is a genuine equidistance locus, and
// it must appear in the edge list -- twice, once per direction, so the two
// cuts' sweep contributions cancel exactly.
TEST(FlatDiagramTest, ACellCanBorderItselfAcrossTheSeam) {
    Diagram diagram(Domain::torus(1.0, 1.0), {Vector2(0.25, 0.5), Vector2(0.75, 0.5)});

    EXPECT_NEAR(diagram.cell(0).area(), 0.5, 1e-9);
    EXPECT_NEAR(diagram.cell(1).area(), 0.5, 1e-9);

    std::set<size_t> neighbors;
    size_t self_edges = 0;

    for (const CellEdgeInfo& edge : diagram.cell_edges(0)) {
        neighbors.insert(edge.neighbor_index);
        self_edges += edge.neighbor_index == 0 ? 1 : 0;
    }

    EXPECT_EQ(neighbors, (std::set<size_t>{0, 1}));
    EXPECT_EQ(self_edges % 2, 0u);
}

TEST(FlatDiagramTest, RebuiltPreservesSitesUpToWrapping) {
    auto diagram = scattered_torus(9, 2.0, 1.0);
    std::vector<Vector3> moved;

    for (size_t index = 0; index < diagram->size(); ++index) {
        moved.push_back(diagram->site_vector(index) + Vector3(2.0, -1.0, 0.0));
    }

    auto rebuilt = diagram->rebuilt(moved);

    ASSERT_EQ(rebuilt->size(), diagram->size());

    for (size_t index = 0; index < diagram->size(); ++index) {
        EXPECT_NEAR((rebuilt->site(index) - diagram->site(index)).norm(), 0.0, 1e-12);
    }
}

// The band of periodic images is laid at a guess from the mean cell width,
// which a clustered scatter defeats: its cells span far more than that, so
// the covering has to be measured and widened before the diagram is right.
TEST(FlatDiagramTest, ClusteredSitesWidenTheBandUntilTheCellsStillPartition) {
    std::vector<Vector2> sites;

    for (const Vector2& site : flat_scatter(64, 1.0, 1.0)) {
        sites.emplace_back(0.02 * site.x(), 0.02 * site.y());
    }

    auto diagram = diagram_of(sites, Domain::torus(1.0, 1.0));

    EXPECT_NEAR(total_cell_area(*diagram), 1.0, 1e-9);
}

// The same grid on the plane: the walls run through the middle of the
// outer cells' would-be neighbours, so every cell is still a 0.5 square.
TEST(PlaneDiagramTest, GridSitesGiveRectangularCellsOfEqualArea) {
    Diagram diagram(Domain::plane(2.0, 1.0), grid_sites());

    for (size_t index = 0; index < diagram.size(); ++index) {
        EXPECT_NEAR(diagram.cell(index).area(), 0.25, 1e-9);
    }
}

TEST(PlaneDiagramTest, CanonicalLeavesPointsWhereTheyAre) {
    Diagram diagram(Domain::plane(2.0, 1.0), {});

    EXPECT_NEAR(diagram.canonical(Vector2(2.5, -0.25)).x(), 2.5, 1e-12);
    EXPECT_NEAR(diagram.canonical(Vector2(2.5, -0.25)).y(), -0.25, 1e-12);
}

TEST(PlaneDiagramTest, CellAreasPartitionTheRectangle) {
    auto diagram = scattered_plane(17, 2.0, 1.0);

    EXPECT_NEAR(total_cell_area(*diagram), 2.0, 1e-9);

    for (size_t index = 0; index < diagram->size(); ++index) {
        for (const Vector2& vertex : diagram->cell(index).vertices()) {
            EXPECT_GE(vertex.x(), -1e-12);
            EXPECT_LE(vertex.x(), 2.0 + 1e-12);
            EXPECT_GE(vertex.y(), -1e-12);
            EXPECT_LE(vertex.y(), 1.0 + 1e-12);
        }
    }
}

// A single site owns the whole rectangle, bounded by walls alone.
TEST(PlaneDiagramTest, ASingleSiteOwnsTheWholeRectangle) {
    Diagram diagram(Domain::plane(2.0, 1.0), {Vector2(0.3, 0.8)});

    EXPECT_NEAR(diagram.cell(0).area(), 2.0, 1e-12);
    EXPECT_TRUE(diagram.cell_edges(0).empty());
}

// Two sites side by side: the bisector runs from the bottom wall to the
// top wall, and both of its endpoint cuts are walls.
TEST(PlaneDiagramTest, ABisectorEndingOnTheWallsIsPinnedByThem) {
    Diagram diagram(Domain::plane(2.0, 1.0), {Vector2(0.5, 0.5), Vector2(1.5, 0.5)});

    std::vector<CellEdgeInfo> edges = diagram.cell_edges(0);
    ASSERT_EQ(edges.size(), 1u);

    const CellEdgeInfo& edge = edges[0];
    EXPECT_EQ(edge.neighbor_index, 1u);
    EXPECT_NEAR(edge.boundary.source().x(), 1.0, 1e-12);
    EXPECT_NEAR(edge.boundary.target().x(), 1.0, 1e-12);
    EXPECT_NEAR(edge.boundary.length(), 1.0, 1e-12);
    EXPECT_TRUE(std::holds_alternative<Wall>(edge.source_cut));
    EXPECT_TRUE(std::holds_alternative<Wall>(edge.target_cut));
    EXPECT_NEAR(diagram.cell(0).area(), 1.0, 1e-12);
    EXPECT_NEAR(diagram.cell(1).area(), 1.0, 1e-12);
}

TEST(PlaneDiagramTest, NeighborPositionsAreTheCanonicalSites) {
    auto diagram = scattered_plane(13, 1.5, 1.0);

    for (size_t index = 0; index < diagram->size(); ++index) {
        for (const CellEdgeInfo& edge : diagram->cell_edges(index)) {
            EXPECT_NE(edge.neighbor_index, index);
            EXPECT_NEAR((edge.neighbor_position - diagram->site(edge.neighbor_index)).norm(), 0.0, 1e-12);
        }
    }
}

// The cylinder wraps its width and walls its height: the grid's cells are
// the same squares, the seam neighbours are period images, and the rim
// cells end at the walls.
TEST(CylinderDiagramTest, GridSitesGiveRectangularCellsOfEqualArea) {
    Diagram diagram(Domain::cylinder(2.0, 1.0), grid_sites());

    for (size_t index = 0; index < diagram.size(); ++index) {
        EXPECT_NEAR(diagram.cell(index).area(), 0.25, 1e-9);
    }
}

TEST(CylinderDiagramTest, WrapsTheWidthAndWallsTheHeight) {
    auto diagram = geometry_art::testing::scattered_diagram(13, Domain::cylinder(1.5, 1.0));
    bool crossed_seam = false;

    EXPECT_NEAR(total_cell_area(*diagram), 1.5, 1e-9);

    for (size_t index = 0; index < diagram->size(); ++index) {
        for (const Vector2& vertex : diagram->cell(index).vertices()) {
            EXPECT_GE(vertex.y(), -1e-12);
            EXPECT_LE(vertex.y(), 1.0 + 1e-12);
        }

        for (const CellEdgeInfo& edge : diagram->cell_edges(index)) {
            Vector2 offset = edge.neighbor_position - diagram->site(edge.neighbor_index);
            crossed_seam = crossed_seam || std::abs(offset.x()) > 1e-9;
            EXPECT_NEAR(offset.y(), 0.0, 1e-12);
        }
    }

    EXPECT_TRUE(crossed_seam);
}

// A site nudged just past a wall keeps a cell: the domain does not move
// with the sites, so the optimizer may step through the wall and back.
TEST(PlaneDiagramTest, ASiteJustOutsideTheWallStillOwnsACell) {
    Diagram diagram(Domain::plane(2.0, 1.0), {Vector2(-0.05, 0.5), Vector2(1.5, 0.5)});

    EXPECT_GT(diagram.cell(0).area(), 0.5);
    EXPECT_NEAR(diagram.cell(0).area() + diagram.cell(1).area(), 2.0, 1e-12);
}
