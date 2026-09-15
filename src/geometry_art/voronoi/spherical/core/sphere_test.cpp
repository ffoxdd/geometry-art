#include <gtest/gtest.h>
#include "sphere.hpp"
#include "../../../types.hpp"
#include "../../../geometry/spherical/helpers.hpp"
#include "../../../testing/assertions/geometric.hpp"

using namespace geometry_art;
using namespace geometry_art::voronoi;
using namespace geometry_art::voronoi::spherical;
using geometry_art::testing::points_approximately_equal;

Sphere create_simple_voronoi_sphere() {
    Sphere sphere;

    sphere.insert(cgal::Point3(1, 0, 0));
    sphere.insert(cgal::Point3(0, 1, 0));
    sphere.insert(cgal::Point3(0, 0, 1));
    sphere.insert(cgal::Point3(-1, 0, 0));

    return sphere;
}

TEST(SphereTest, DefaultConstructorCreatesEmptySphere) {
    Sphere sphere;
    EXPECT_EQ(sphere.size(), 0);
}

TEST(SphereTest, InsertIncreasesSize) {
    Sphere sphere;
    EXPECT_EQ(sphere.size(), 0);

    sphere.insert(cgal::Point3(1, 0, 0));
    EXPECT_EQ(sphere.size(), 1);

    sphere.insert(cgal::Point3(0, 1, 0));
    EXPECT_EQ(sphere.size(), 2);

    sphere.insert(cgal::Point3(0, 0, 1));
    EXPECT_EQ(sphere.size(), 3);
}

TEST(SphereTest, SiteReturnsInsertedPoint) {
    Sphere sphere;
    cgal::Point3 point(1, 0, 0);
    sphere.insert(point);

    cgal::Point3 retrieved = sphere.site(0);
    EXPECT_TRUE(points_approximately_equal(point, retrieved));
}

TEST(SphereTest, SiteReturnsCorrectPointsForMultipleSites) {
    Sphere sphere;
    cgal::Point3 p1(1, 0, 0);
    cgal::Point3 p2(0, 1, 0);
    cgal::Point3 p3(0, 0, 1);

    sphere.insert(p1);
    sphere.insert(p2);
    sphere.insert(p3);

    EXPECT_TRUE(points_approximately_equal(sphere.site(0), p1));
    EXPECT_TRUE(points_approximately_equal(sphere.site(1), p2));
    EXPECT_TRUE(points_approximately_equal(sphere.site(2), p3));
}

TEST(SphereTest, UpdateSiteChangesPosition) {
    Sphere sphere;
    cgal::Point3 original(1, 0, 0);
    cgal::Point3 updated(0, 1, 0);

    sphere.insert(original);
    EXPECT_TRUE(points_approximately_equal(sphere.site(0), original));

    sphere.update_site(0, updated);
    EXPECT_TRUE(points_approximately_equal(sphere.site(0), updated));
}

TEST(SphereTest, CellArcsReturnsNonEmptyRange) {
    Sphere sphere = create_simple_voronoi_sphere();

    size_t arc_count = 0;
    for (const auto &cell : sphere.cells()) {
        arc_count += cell.arcs().size();
    }

    EXPECT_GT(arc_count, 0);
}

TEST(SphereTest, CellArcsIsEmptyForEmptySphere) {
    Sphere sphere;

    size_t arc_count = 0;
    for (const auto &cell : sphere.cells()) {
        arc_count += cell.arcs().size();
    }

    EXPECT_EQ(arc_count, 0);
}

TEST(SphereTest, CellArcsIsEmptyForSinglePoint) {
    Sphere sphere;
    sphere.insert(cgal::Point3(1, 0, 0));

    size_t arc_count = 0;
    for (const auto &cell : sphere.cells()) {
        arc_count += cell.arcs().size();
    }

    EXPECT_EQ(arc_count, 0);
}

TEST(SphereTest, DualCellsReturnsPolygons) {
    Sphere sphere = create_simple_voronoi_sphere();

    size_t cell_count = 0;
    for (const auto &cell : sphere.cells()) {
        // Constructor already validates with CGAL_precondition
        cell_count++;
    }

    EXPECT_EQ(cell_count, sphere.size());
}

TEST(SphereTest, DualCellsIsEmptyForEmptySphere) {
    Sphere sphere;

    size_t cell_count = 0;
    for (const auto &cell : sphere.cells()) {
        cell_count++;
    }

    EXPECT_EQ(cell_count, 0);
}

TEST(SphereTest, MultipleUpdatesToSameSite) {
    Sphere sphere;
    sphere.insert(cgal::Point3(1, 0, 0));

    cgal::Point3 pos1 = cgal::to_point(VectorS2(0, 1, 0).normalized());
    cgal::Point3 pos2 = cgal::to_point(VectorS2(0, 0, 1).normalized());
    cgal::Point3 pos3 = cgal::to_point(VectorS2(1, 1, 1).normalized());

    sphere.update_site(0, pos1);
    EXPECT_TRUE(points_approximately_equal(sphere.site(0), pos1));

    sphere.update_site(0, pos2);
    EXPECT_TRUE(points_approximately_equal(sphere.site(0), pos2));

    sphere.update_site(0, pos3);
    EXPECT_TRUE(points_approximately_equal(sphere.site(0), pos3));
}

TEST(SphereTest, UpdateSitePreservesOtherSites) {
    Sphere sphere = create_simple_voronoi_sphere();
    cgal::Point3 p0 = sphere.site(0);
    cgal::Point3 p1 = sphere.site(1);
    cgal::Point3 p2 = sphere.site(2);
    cgal::Point3 p3 = sphere.site(3);

    // Update one site to a different but valid position
    cgal::Point3 updated_p1 = cgal::to_point(VectorS2(0.707, 0.707, 0).normalized());
    sphere.update_site(1, updated_p1);

    // Other sites should remain unchanged
    EXPECT_TRUE(points_approximately_equal(sphere.site(0), p0));
    EXPECT_TRUE(points_approximately_equal(sphere.site(1), updated_p1));
    EXPECT_TRUE(points_approximately_equal(sphere.site(2), p2));
    EXPECT_TRUE(points_approximately_equal(sphere.site(3), p3));
}

TEST(SphereTest, ArcsReturnsAllCellArcs) {
    Sphere sphere = create_simple_voronoi_sphere();

    size_t arcs_via_cells = 0;
    for (const auto &cell : sphere.cells()) {
        arcs_via_cells += cell.arcs().size();
    }

    size_t arcs_via_arcs = 0;
    for (const auto &arc : sphere.arcs()) {
        (void)arc;
        arcs_via_arcs++;
    }

    EXPECT_EQ(arcs_via_arcs, arcs_via_cells);
    EXPECT_GT(arcs_via_arcs, 0);
}

TEST(SphereTest, ArcsIsEmptyForEmptySphere) {
    Sphere sphere;

    size_t arc_count = 0;
    for (const auto &arc : sphere.arcs()) {
        (void)arc;
        arc_count++;
    }

    EXPECT_EQ(arc_count, 0);
}

TEST(SphereTest, ArcsIsEmptyForSinglePoint) {
    Sphere sphere;
    sphere.insert(cgal::Point3(1, 0, 0));

    size_t arc_count = 0;
    for (const auto &arc : sphere.arcs()) {
        (void)arc;
        arc_count++;
    }

    EXPECT_EQ(arc_count, 0);
}


TEST(SphereTest, CellEdgeArcsLieOnTheBisectorsOfTheirSites) {
    Sphere sphere;
    sphere.insert(cgal::to_point(VectorS2(1, 0.1, 0.2).normalized()));
    sphere.insert(cgal::to_point(VectorS2(-0.3, 1, 0.1).normalized()));
    sphere.insert(cgal::to_point(VectorS2(0.2, -0.4, 1).normalized()));
    sphere.insert(cgal::to_point(VectorS2(-1, -0.6, -0.5).normalized()));
    sphere.insert(cgal::to_point(VectorS2(0.5, -1, -0.3).normalized()));

    for (size_t k = 0; k < sphere.size(); ++k) {
        Vector3 site = to_vector3(sphere.site(k));

        for (const auto& edge : sphere.cell_edges(k)) {
            Vector3 neighbor = to_vector3(sphere.site(edge.neighbor_index));
            VectorS2 bisector = VectorS2(site - neighbor).normalized();
            EXPECT_NEAR((edge.arc.normal() - bisector).norm(), 0.0, 1e-12);
            EXPECT_LT(edge.arc.length(), M_PI);
        }
    }
}

TEST(SphereTest, CellArcsRunCounterClockwiseAboutTheirCell) {
    Sphere sphere = create_simple_voronoi_sphere();

    for (size_t index = 0; index < sphere.size(); ++index) {
        VectorS2 own = to_vector_s2(sphere.site(index));

        for (const CellEdgeInfo& edge : sphere.cell_edges(index)) {
            VectorS2 neighbor = to_vector_s2(sphere.site(edge.neighbor_index));
            EXPECT_GT(edge.arc.normal().dot(own - neighbor), 0.0) << "cell " << index;
        }
    }
}

TEST(SphereTest, OffsetCellStaysTheInsetDistanceInsideEveryBisector) {
    Sphere sphere = create_simple_voronoi_sphere();
    double inset = 0.1;

    for (size_t index = 0; index < sphere.size(); ++index) {
        CapPolygon region = sphere.offset_cell(index, inset);
        VectorS2 own = to_vector_s2(sphere.site(index));

        ASSERT_GE(region.size(), 2u) << "cell " << index;

        for (const CapPolygon::Edge& edge : region.edges()) {
            EXPECT_NEAR(edge.cap.offset, std::sin(inset), 1e-12);

            for (const CellEdgeInfo& cell_edge : sphere.cell_edges(index)) {
                VectorS2 neighbor = to_vector_s2(sphere.site(cell_edge.neighbor_index));
                EXPECT_GE((own - neighbor).normalized().dot(edge.source), std::sin(inset) - 1e-9) << "cell " << index;
            }
        }
    }
}

TEST(SphereTest, OffsetCellOfASingleSiteIsEmpty) {
    Sphere sphere;
    sphere.insert(cgal::Point3(1, 0, 0));

    EXPECT_TRUE(sphere.offset_cell(0, 0.1).empty());
}
