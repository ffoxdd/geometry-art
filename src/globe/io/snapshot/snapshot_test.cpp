#include "snapshot.hpp"
#include "../../fields/spherical/polynomial_field.hpp"
#include "../../testing/explicit_diagram.hpp"
#include <gtest/gtest.h>
#include <cmath>

using namespace globe;
using fields::spherical::PolynomialField;
using io::snapshot::Snapshot;
using io::snapshot::capture;
using globe::testing::ExplicitDiagram;

TEST(SnapshotTest, CapturesEveryCellOfTheDiagram) {
    Snapshot snapshot = capture(ExplicitDiagram::hemispheres(), PolynomialField::constant(1.0));

    ASSERT_EQ(snapshot.cells.size(), 2u);
    EXPECT_EQ(snapshot.geometry, "sphere");
    EXPECT_EQ(snapshot.cells[0].site_index, 0u);
    EXPECT_EQ(snapshot.cells[1].site_index, 1u);
}

TEST(SnapshotTest, HemispheresCarryHalfTheMassEach) {
    Snapshot snapshot = capture(ExplicitDiagram::hemispheres(), PolynomialField::constant(1.0));

    EXPECT_NEAR(snapshot.total_mass, 4.0 * M_PI, 1e-12);
    EXPECT_NEAR(snapshot.cells[0].mass, 2.0 * M_PI, 1e-12);
    EXPECT_NEAR(snapshot.cells[1].mass, 2.0 * M_PI, 1e-12);
    EXPECT_NEAR(snapshot.cells[0].area, 2.0 * M_PI, 1e-12);
}

TEST(SnapshotTest, EqualCapacitiesLeaveNoRelativeError) {
    Snapshot snapshot = capture(ExplicitDiagram::hemispheres(), PolynomialField::constant(1.0));

    EXPECT_NEAR(snapshot.target_mass(), 2.0 * M_PI, 1e-12);
    EXPECT_NEAR(snapshot.relative_rms_capacity_error(), 0.0, 1e-12);
}

TEST(SnapshotTest, UnequalCapacitiesShowAsRelativeError) {
    // A density leaning north puts more mass in the northern cell.
    Snapshot snapshot = capture(
        ExplicitDiagram::hemispheres(),
        PolynomialField::linear(2.0, Vector3(0.0, 0.0, 1.0))
    );

    EXPECT_GT(snapshot.cells[0].mass, snapshot.cells[1].mass);
    EXPECT_GT(snapshot.relative_rms_capacity_error(), 0.0);
    EXPECT_NEAR(snapshot.cells[0].mass + snapshot.cells[1].mass, snapshot.total_mass, 1e-12);
}

TEST(SnapshotTest, BoundaryAndNeighborsComeFromTheDiagram) {
    Snapshot snapshot = capture(ExplicitDiagram::hemispheres(), PolynomialField::constant(1.0));

    EXPECT_EQ(snapshot.cells[0].boundary.size(), 4u);
    EXPECT_EQ(snapshot.cells[0].neighbors, std::vector<size_t>({1, 1, 1, 1}));

    for (const Vector3& point : snapshot.cells[0].boundary) {
        EXPECT_NEAR(point.z(), 0.0, 1e-12);
        EXPECT_NEAR(point.norm(), 1.0, 1e-12);
    }
}

TEST(SnapshotTest, ReportsNothingForAnEmptyDiagram) {
    Snapshot snapshot;

    EXPECT_NEAR(snapshot.target_mass(), 0.0, 1e-15);
    EXPECT_NEAR(snapshot.relative_rms_capacity_error(), 0.0, 1e-15);
}
