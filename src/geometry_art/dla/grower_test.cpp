#include "gtest/gtest.h"
#include "grower.hpp"
#include "aggregate.hpp"
#include "callback.hpp"
#include "parameters.hpp"
#include "particle.hpp"
#include "walker.hpp"
#include "direction_sampler/uniform_direction_sampler.hpp"
#include "spawn_shell/harmonic_spawn_shell.hpp"
#include "spawn_shell/killing_spawn_shell.hpp"
#include "../math/interval_sampler/uniform_interval_sampler.hpp"
#include "../testing/macros.hpp"
#include "../testing/mocks/direction_sampler.hpp"
#include "../testing/mocks/spawn_shell.hpp"
#include "../types.hpp"
#include <Eigen/Eigenvalues>
#include <cmath>
#include <cstddef>
#include <optional>
#include <vector>

using namespace geometry_art::dla;
using geometry_art::UniformIntervalSampler;
using geometry_art::Vector3;
using geometry_art::testing::mocks::MockDirectionSampler;
using geometry_art::testing::mocks::MockSpawnShell;

namespace {

Parameters chain_parameters() {
    Parameters parameters;
    parameters.particle_count = 3;
    parameters.particle_radius = 1.0;
    parameters.overlap = 0.01;
    return parameters;
}

Aggregate<> seeded_aggregate() {
    return Aggregate<>(1.0, Particle{Vector3::Zero(), std::nullopt});
}

} // namespace

TEST(GrowerTest, GrowsAChainToTheRequestedCount) {
    Parameters parameters = chain_parameters();

    Walker<MockDirectionSampler, MockSpawnShell> walker(
        parameters,
        MockDirectionSampler({Vector3(0.0, 0.0, -1.0)}),
        MockSpawnShell(Vector3(0.0, 0.0, 5.0))
    );

    size_t reports = 0;
    Callback<> callback = [&reports](const Aggregate<>&) { ++reports; };

    Grower<RTreeParticleIndex, MockDirectionSampler, MockSpawnShell> grower(
        parameters,
        seeded_aggregate(),
        std::move(walker),
        callback
    );

    Aggregate<> aggregate = grower.grow();

    ASSERT_EQ(aggregate.size(), 3u);
    EXPECT_EQ(reports, 2u);

    EXPECT_NEAR((aggregate.particles()[1].center - Vector3(0.0, 0.0, 1.99)).norm(), 0.0, 1e-9);
    EXPECT_NEAR((aggregate.particles()[2].center - Vector3(0.0, 0.0, 3.98)).norm(), 0.0, 1e-9);
    EXPECT_EQ(aggregate.particles()[1].parent, std::optional<std::size_t>(0));
    EXPECT_EQ(aggregate.particles()[2].parent, std::optional<std::size_t>(1));
}

namespace {

Parameters growth_parameters(size_t particle_count) {
    Parameters parameters;
    parameters.particle_count = particle_count;
    parameters.particle_radius = 1.0;
    parameters.overlap = 0.01;
    parameters.spawn_margin = 2.0;
    return parameters;
}

Aggregate<> grow_harmonic(const Parameters& parameters, unsigned int seed) {
    Walker<UniformDirectionSampler<>, HarmonicSpawnShell<>> walker(
        parameters,
        UniformDirectionSampler<>(UniformIntervalSampler(seed)),
        HarmonicSpawnShell<>(
            UniformDirectionSampler<>(UniformIntervalSampler(seed + 1)),
            UniformIntervalSampler(seed + 2)
        )
    );

    Grower<> grower(parameters, seeded_aggregate(), std::move(walker), noop_callback());
    return grower.grow();
}

Aggregate<> grow_killing(const Parameters& parameters, unsigned int seed, double kill_factor) {
    Walker<UniformDirectionSampler<>, KillingSpawnShell<>> walker(
        parameters,
        UniformDirectionSampler<>(UniformIntervalSampler(seed)),
        KillingSpawnShell<>(kill_factor, UniformDirectionSampler<>(UniformIntervalSampler(seed + 1)))
    );

    Grower<RTreeParticleIndex, UniformDirectionSampler<>, KillingSpawnShell<>> grower(
        parameters,
        seeded_aggregate(),
        std::move(walker),
        noop_callback()
    );

    return grower.grow();
}

double mean_center_norm(const Aggregate<>& aggregate) {
    double total = 0.0;

    for (const Particle& particle : aggregate.particles()) {
        total += particle.center.norm();
    }

    return total / static_cast<double>(aggregate.size());
}

} // namespace

TEST(GrowerTest, EXPENSIVE_GrownAggregateIsConnectedAndBarelyOverlapping) {
    REQUIRE_EXPENSIVE();

    Parameters parameters = growth_parameters(400);
    Aggregate<> aggregate = grow_harmonic(parameters, 7);

    ASSERT_EQ(aggregate.size(), 400u);
    EXPECT_GT(aggregate.reach(), 2.0);

    const std::vector<Particle>& particles = aggregate.particles();
    double stuck_distance = aggregate.contact_distance() - parameters.overlap_distance();

    for (size_t index = 1; index < particles.size(); ++index) {
        ASSERT_TRUE(particles[index].parent.has_value());
        ASSERT_LT(*particles[index].parent, index);

        double parent_distance = (particles[index].center - particles[*particles[index].parent].center).norm();
        EXPECT_NEAR(parent_distance, stuck_distance, 1e-9);
    }

    double interpenetration_floor = aggregate.contact_distance() - 2.0 * parameters.overlap_distance() - 1e-9;

    for (size_t a = 0; a < particles.size(); ++a) {
        for (size_t b = a + 1; b < particles.size(); ++b) {
            EXPECT_GE((particles[a].center - particles[b].center).norm(), interpenetration_floor);
        }
    }
}

// The reference does nothing clever: a roomy spawn shell and a kill radius
// so far out that killing early walkers cannot matter. The production
// configuration -- tight shell, exact harmonic returns -- must grow the same
// structure.
TEST(GrowerTest, EXPENSIVE_TightHarmonicShellMatchesAConservativeReference) {
    REQUIRE_EXPENSIVE();

    Parameters tight = growth_parameters(300);

    Parameters roomy = growth_parameters(300);
    roomy.spawn_margin = 10.0;

    double harmonic_mean = mean_center_norm(grow_harmonic(tight, 11));
    double reference_mean = mean_center_norm(grow_killing(roomy, 23, 30.0));

    EXPECT_NEAR(harmonic_mean / reference_mean, 1.0, 0.25);
}

// A tripwire for directional bias: shape fluctuations are large in DLA, but
// a lopsided sampler or a skewed return kernel would push the center of mass
// off the seed and stretch the principal axes far beyond these bounds.
TEST(GrowerTest, EXPENSIVE_GrowthIsStatisticallyIsotropic) {
    REQUIRE_EXPENSIVE();

    for (unsigned int seed : {3u, 5u, 7u}) {
        Aggregate<> aggregate = grow_harmonic(growth_parameters(600), seed);

        Vector3 center_of_mass = Vector3::Zero();

        for (const Particle& particle : aggregate.particles()) {
            center_of_mass += particle.center;
        }

        center_of_mass /= static_cast<double>(aggregate.size());

        geometry_art::Matrix3 covariance = geometry_art::Matrix3::Zero();

        for (const Particle& particle : aggregate.particles()) {
            Vector3 offset = particle.center - center_of_mass;
            covariance += offset * offset.transpose();
        }

        covariance /= static_cast<double>(aggregate.size());

        Eigen::SelfAdjointEigenSolver<geometry_art::Matrix3> solver(covariance);
        Vector3 axis_lengths = solver.eigenvalues().cwiseSqrt();
        double gyration_radius = std::sqrt(solver.eigenvalues().sum());

        EXPECT_LT(center_of_mass.norm() / gyration_radius, 0.5) << "seed " << seed;
        EXPECT_LT(axis_lengths(2) / axis_lengths(0), 2.2) << "seed " << seed;
    }
}
