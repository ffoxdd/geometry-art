#include "gtest/gtest.h"
#include "walker.hpp"
#include "aggregate.hpp"
#include "parameters.hpp"
#include "particle.hpp"
#include "../testing/mocks/direction_sampler.hpp"
#include "../testing/mocks/spawn_shell.hpp"
#include "../types.hpp"
#include <optional>
#include <vector>

using namespace geometry_art::dla;
using geometry_art::Vector3;
using geometry_art::testing::mocks::MockDirectionSampler;
using geometry_art::testing::mocks::MockSpawnShell;

using MockedWalker = Walker<MockDirectionSampler, MockSpawnShell>;

namespace {

Parameters test_parameters() {
    Parameters parameters;
    parameters.particle_radius = 1.0;
    parameters.overlap = 0.01;
    return parameters;
}

Aggregate<> seeded_aggregate() {
    return Aggregate<>(1.0, Particle{Vector3::Zero(), std::nullopt});
}

} // namespace

TEST(WalkerTest, WalksStraightDownOntoTheSeed) {
    Aggregate<> aggregate = seeded_aggregate();

    MockedWalker walker(
        test_parameters(),
        MockDirectionSampler({Vector3(0.0, 0.0, -1.0)}),
        MockSpawnShell(Vector3(0.0, 0.0, 5.0))
    );

    Particle particle = walker.settle(aggregate);

    EXPECT_NEAR((particle.center - Vector3(0.0, 0.0, 1.99)).norm(), 0.0, 1e-9);
    EXPECT_EQ(particle.parent, std::optional<std::size_t>(0));
}

TEST(WalkerTest, SticksToTheParticleItHits) {
    Aggregate<> aggregate = seeded_aggregate();
    aggregate.freeze(Particle{Vector3(0.0, 0.0, 1.99), 0});

    MockedWalker walker(
        test_parameters(),
        MockDirectionSampler({Vector3(0.0, 0.0, -1.0)}),
        MockSpawnShell(Vector3(0.0, 0.0, 6.0))
    );

    Particle particle = walker.settle(aggregate);

    EXPECT_NEAR((particle.center - Vector3(0.0, 0.0, 3.98)).norm(), 0.0, 1e-9);
    EXPECT_EQ(particle.parent, std::optional<std::size_t>(1));
}
