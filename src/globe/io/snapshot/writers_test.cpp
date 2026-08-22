#include "json_writer.hpp"
#include "svg_writer.hpp"
#include "../../fields/spherical/polynomial_field.hpp"
#include "../../generators/spherical/fibonacci_point_generator.hpp"
#include "../../testing/explicit_diagram.hpp"
#include "../../testing/golden_file.hpp"
#include "../../voronoi/spherical/core/random_builder.hpp"
#include <gtest/gtest.h>
#include <memory>

using namespace globe;
using fields::spherical::PolynomialField;
using generators::spherical::FibonacciPointGenerator;
using io::snapshot::JsonWriter;
using io::snapshot::Snapshot;
using io::snapshot::SvgWriter;
using io::snapshot::capture;
using voronoi::spherical::RandomBuilder;
using globe::testing::ExplicitDiagram;
using globe::testing::GoldenFile;

namespace {

// Fibonacci placement is deterministic, so the snapshot is fixed without an
// iterative solve whose last digits could differ between platforms.
Snapshot fibonacci_snapshot(int count) {
    auto sphere = RandomBuilder<FibonacciPointGenerator>(FibonacciPointGenerator()).build(count);
    return capture(*sphere, PolynomialField::constant(1.0));
}

std::string golden(const std::string& name) {
    return "src/globe/io/snapshot/testdata/" + name;
}

} // namespace

TEST(JsonWriterTest, WritesTheHemispheresSnapshot) {
    Snapshot snapshot = capture(ExplicitDiagram::hemispheres(), PolynomialField::constant(1.0));

    GoldenFile(golden("hemispheres.json")).expect_matches(JsonWriter().to_string(snapshot));
}

TEST(JsonWriterTest, WritesAFibonacciSnapshot) {
    GoldenFile(golden("fibonacci_24.json")).expect_matches(JsonWriter().to_string(fibonacci_snapshot(24)));
}

TEST(JsonWriterTest, RoundsNegativeZeroToZero) {
    Snapshot snapshot;
    snapshot.geometry = "sphere";
    snapshot.cells.push_back(Snapshot::Cell{0, Vector3(-0.0, 0.0, 1.0), {}, {}, 0.0, 0.0});

    EXPECT_EQ(JsonWriter().to_string(snapshot).find("-0"), std::string::npos);
}

TEST(SvgWriterTest, DrawsAFibonacciSnapshot) {
    GoldenFile(golden("fibonacci_24.svg")).expect_matches(SvgWriter().to_string(fibonacci_snapshot(24)));
}

TEST(SvgWriterTest, DrawsOnlyTheFacingHemisphere) {
    Snapshot snapshot = fibonacci_snapshot(24);
    std::string drawing = SvgWriter().to_string(snapshot);

    size_t drawn = 0;
    size_t position = drawing.find("<path");

    while (position != std::string::npos) {
        ++drawn;
        position = drawing.find("<path", position + 1);
    }

    size_t facing = 0;

    for (const Snapshot::Cell& cell : snapshot.cells) {
        facing += cell.site.z() > 0.0 ? 1 : 0;
    }

    EXPECT_EQ(drawn, facing);
    EXPECT_LT(drawn, snapshot.cells.size());
}
