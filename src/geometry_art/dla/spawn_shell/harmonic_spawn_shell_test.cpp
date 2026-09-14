#include "gtest/gtest.h"
#include "spawn_shell.hpp"
#include "harmonic_spawn_shell.hpp"
#include "../../testing/mocks/direction_sampler.hpp"
#include "../../testing/mocks/interval_sampler.hpp"
#include "../../types.hpp"
#include <vector>

using namespace geometry_art::dla;
using geometry_art::Vector3;
using geometry_art::testing::mocks::MockDirectionSampler;
using geometry_art::testing::mocks::MockIntervalSampler;

using MockedShell = HarmonicSpawnShell<MockDirectionSampler, MockIntervalSampler>;

TEST(HarmonicSpawnShellTest, SpawnScalesTheDirectionToTheShell) {
    MockedShell shell(MockDirectionSampler({Vector3(1.0, 0.0, 0.0)}), MockIntervalSampler({0.0}));

    Vector3 position = shell.spawn(3.0);

    EXPECT_NEAR((position - Vector3(3.0, 0.0, 0.0)).norm(), 0.0, 1e-12);
}

TEST(HarmonicSpawnShellTest, ConfineLeavesTheInteriorAlone) {
    MockedShell shell(MockDirectionSampler({Vector3(1.0, 0.0, 0.0)}), MockIntervalSampler({0.0}));
    Vector3 position(0.0, 0.0, 1.0);

    EXPECT_EQ(shell.confine(position, 2.0), position);
}

TEST(HarmonicSpawnShellTest, ConfineLandsAReturningWalkerByTheHarmonicKernel) {
    MockedShell shell(
        MockDirectionSampler({Vector3(1.0, 0.0, 0.0)}),
        MockIntervalSampler({0.0, 1.0, 0.0})
    );

    Vector3 landing = shell.confine(Vector3(0.0, 0.0, 4.0), 2.0);

    EXPECT_NEAR((landing - Vector3(0.0, 0.0, 2.0)).norm(), 0.0, 1e-9);
}

TEST(HarmonicSpawnShellTest, ConfineRespawnsAnEscapedWalker) {
    MockedShell shell(
        MockDirectionSampler({Vector3(0.0, 1.0, 0.0)}),
        MockIntervalSampler({0.9})
    );

    Vector3 respawned = shell.confine(Vector3(0.0, 0.0, 4.0), 2.0);

    EXPECT_NEAR((respawned - Vector3(0.0, 2.0, 0.0)).norm(), 0.0, 1e-12);
}
