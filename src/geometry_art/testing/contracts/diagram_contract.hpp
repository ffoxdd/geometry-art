#ifndef GEOMETRY_ART_TESTING_CONTRACTS_DIAGRAM_CONTRACT_HPP_
#define GEOMETRY_ART_TESTING_CONTRACTS_DIAGRAM_CONTRACT_HPP_

#include "../../types.hpp"
#include <gtest/gtest.h>
#include <cstddef>
#include <memory>
#include <vector>

namespace geometry_art::testing::contracts {

struct EdgeRecord {
    size_t neighbor_index;
    Vector3 neighbor_position;
    Vector3 source;
    Vector3 target;
};

// Laws any Voronoi diagram must satisfy, stated against a traits type:
//
//   using DiagramType = ...;
//   static std::unique_ptr<DiagramType> scattered(size_t count);
//   static double domain_area();
//   static double cell_area(const DiagramType&, size_t index);
//   static Vector3 site(const DiagramType&, size_t index);
//   static std::vector<EdgeRecord> edges(const DiagramType&, size_t index);
//
// Each law reads one edge snapshot per cell -- a diagram's edge ordering
// is not part of the contract, so no law may correlate separate queries.
// Positions are reported in the cell's own chart, so equidistance is
// checked with the ambient metric on both geometries. Example-based tests
// stay with each implementation; these are only the laws.
template<typename Traits>
class DiagramContract : public ::testing::Test {};

TYPED_TEST_SUITE_P(DiagramContract);

TYPED_TEST_P(DiagramContract, CellsPartitionTheDomain) {
    auto diagram = TypeParam::scattered(17);
    double total = 0.0;

    for (size_t index = 0; index < diagram->size(); ++index) {
        double area = TypeParam::cell_area(*diagram, index);
        EXPECT_GT(area, 0.0);
        total += area;
    }

    EXPECT_NEAR(total, TypeParam::domain_area(), 1e-8 * TypeParam::domain_area());
}

TYPED_TEST_P(DiagramContract, NeighborsAreSymmetric) {
    auto diagram = TypeParam::scattered(13);

    for (size_t index = 0; index < diagram->size(); ++index) {
        for (const EdgeRecord& edge : TypeParam::edges(*diagram, index)) {
            ASSERT_LT(edge.neighbor_index, diagram->size());
            bool found = false;

            for (const EdgeRecord& back : TypeParam::edges(*diagram, edge.neighbor_index)) {
                found = found || back.neighbor_index == index;
            }

            EXPECT_TRUE(found) << index << " -> " << edge.neighbor_index;
        }
    }
}

TYPED_TEST_P(DiagramContract, BoundariesAreEquidistantFromBothSites) {
    auto diagram = TypeParam::scattered(11);

    for (size_t index = 0; index < diagram->size(); ++index) {
        Vector3 own = TypeParam::site(*diagram, index);

        for (const EdgeRecord& edge : TypeParam::edges(*diagram, index)) {
            for (const Vector3& endpoint : {edge.source, edge.target}) {
                EXPECT_NEAR((endpoint - own).norm(), (endpoint - edge.neighbor_position).norm(), 1e-8);
            }
        }
    }
}

TYPED_TEST_P(DiagramContract, SelfEdgesComeInPairs) {
    auto diagram = TypeParam::scattered(9);

    for (size_t index = 0; index < diagram->size(); ++index) {
        size_t self_edges = 0;

        for (const EdgeRecord& edge : TypeParam::edges(*diagram, index)) {
            self_edges += edge.neighbor_index == index ? 1 : 0;
        }

        EXPECT_EQ(self_edges % 2, 0u) << "cell " << index;
    }
}

REGISTER_TYPED_TEST_SUITE_P(
    DiagramContract,
    CellsPartitionTheDomain,
    NeighborsAreSymmetric,
    BoundariesAreEquidistantFromBothSites,
    SelfEdgesComeInPairs
);

} // namespace geometry_art::testing::contracts

#endif //GEOMETRY_ART_TESTING_CONTRACTS_DIAGRAM_CONTRACT_HPP_
