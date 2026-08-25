#ifndef GEOMETRY_ART_TESTING_CONTRACTS_STATE_CONTRACT_HPP_
#define GEOMETRY_ART_TESTING_CONTRACTS_STATE_CONTRACT_HPP_

#include "../../types.hpp"
#include "../../voronoi/state.hpp"
#include <gtest/gtest.h>
#include <cstddef>
#include <memory>

namespace geometry_art::testing::contracts {

// Laws the assembled diagram state must satisfy for any geometry and any
// admissible field, stated against a traits type:
//
//   static voronoi::DiagramState scattered_state(size_t count);
//   static double total_mass();
//
// The relative-moment identity is the sharp one: the moments about the two
// sites of one bisector differ by exactly the separation times the mass,
// because (x - s_own) - (x - s_neighbor) is the constant site difference.
template<typename Traits>
class StateContract : public ::testing::Test {};

TYPED_TEST_SUITE_P(StateContract);

TYPED_TEST_P(StateContract, CapacitiesPartitionTheTotalMass) {
    voronoi::DiagramState state = TypeParam::scattered_state(14);
    double total = 0.0;

    for (const voronoi::CellState& cell : state.cells) {
        EXPECT_GT(cell.mass, 0.0);
        EXPECT_GE(cell.squared_norm_moment, 0.0);
        total += cell.mass;
    }

    EXPECT_NEAR(total, TypeParam::total_mass(), 1e-8 * TypeParam::total_mass());
}

TYPED_TEST_P(StateContract, RelativeMomentsSatisfyTheSeparationIdentity) {
    voronoi::DiagramState state = TypeParam::scattered_state(12);

    for (size_t k = 0; k < state.edges.size(); ++k) {
        for (const voronoi::EdgeState& edge : state.edges[k]) {
            EXPECT_GE(edge.mass, 0.0);
            EXPECT_GT(edge.separation, 0.0);

            Vector3 difference = edge.moment_about_own - edge.moment_about_neighbor;

            EXPECT_NEAR(difference.norm(), edge.separation * edge.mass, 1e-9 * (1.0 + edge.mass));
        }
    }
}

REGISTER_TYPED_TEST_SUITE_P(
    StateContract,
    CapacitiesPartitionTheTotalMass,
    RelativeMomentsSatisfyTheSeparationIdentity
);

} // namespace geometry_art::testing::contracts

#endif //GEOMETRY_ART_TESTING_CONTRACTS_STATE_CONTRACT_HPP_
