#include "gtest/gtest.h"
#include "spawn_shell.hpp"
#include "killing_spawn_shell.hpp"
#include "../../testing/mocks/direction_sampler.hpp"
#include "../../types.hpp"
#include <vector>

using namespace geometry_art::dla;
using geometry_art::Vector3;
using geometry_art::testing::mocks::MockDirectionSampler;

using MockedShell = KillingSpawnShell<MockDirectionSampler>;

TEST(KillingSpawnShellTest, ConfineLeavesTheInteriorAlone) {
    MockedShell shell(3.0, MockDirectionSampler({Vector3(1.0, 0.0, 0.0)}));
    Vector3 position(0.0, 0.0, 1.0);

    EXPECT_EQ(shell.confine(position, 2.0), position);
}

TEST(KillingSpawnShellTest, ConfineLetsAWandererRoamInsideTheKillRadius) {
    MockedShell shell(3.0, MockDirectionSampler({Vector3(1.0, 0.0, 0.0)}));
    Vector3 position(0.0, 0.0, 5.0);

    EXPECT_EQ(shell.confine(position, 2.0), position);
}

TEST(KillingSpawnShellTest, ConfineRespawnsBeyondTheKillRadius) {
    MockedShell shell(3.0, MockDirectionSampler({Vector3(1.0, 0.0, 0.0)}));

    Vector3 respawned = shell.confine(Vector3(0.0, 0.0, 7.0), 2.0);

    EXPECT_NEAR((respawned - Vector3(2.0, 0.0, 0.0)).norm(), 0.0, 1e-12);
}
