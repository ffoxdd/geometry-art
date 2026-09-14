#include "gtest/gtest.h"
#include "aggregate_snapshot.hpp"
#include "aggregate_json_writer.hpp"
#include "aggregate_svg_writer.hpp"
#include "../../testing/golden_file.hpp"
#include "../../types.hpp"
#include <optional>
#include <string>

using namespace geometry_art::io::snapshot;
using geometry_art::Vector3;
using geometry_art::testing::GoldenFile;

namespace {

std::string golden(const std::string& name) {
    return "src/geometry_art/io/snapshot/testdata/" + name;
}

AggregateSnapshot chain_snapshot() {
    AggregateSnapshot snapshot;
    snapshot.particle_radius = 1.0;
    snapshot.reach = 5.98;

    snapshot.particles = {
        {Vector3(0.0, 0.0, 0.0), std::nullopt},
        {Vector3(0.0, 0.0, 1.99), 0},
        {Vector3(0.5, 0.0, 3.98), 1},
    };

    return snapshot;
}

} // namespace

TEST(AggregateJsonWriterTest, SeedParticleHasANullParent) {
    std::string json = AggregateJsonWriter().to_string(chain_snapshot());

    EXPECT_NE(json.find("\"parent\": null"), std::string::npos);
}

TEST(AggregateJsonWriterTest, MatchesGoldenFile) {
    GoldenFile(golden("aggregate_chain.json")).expect_matches(
        AggregateJsonWriter().to_string(chain_snapshot())
    );
}

TEST(AggregateSvgWriterTest, MatchesGoldenFile) {
    GoldenFile(golden("aggregate_chain.svg")).expect_matches(
        AggregateSvgWriter().to_string(chain_snapshot())
    );
}
