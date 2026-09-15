#include "json_reader.hpp"
#include "json_writer.hpp"
#include "snapshot.hpp"
#include "../../types.hpp"
#include <gtest/gtest.h>
#include <stdexcept>
#include <string>

namespace geometry_art::io::snapshot {
namespace {

std::string testdata(const std::string& name) {
    return std::string(GEOMETRY_ART_SOURCE_DIR) + "/src/geometry_art/io/snapshot/testdata/" + name;
}

Snapshot two_cells() {
    Snapshot snapshot;
    snapshot.geometry = "plane";
    snapshot.width = 2.0;
    snapshot.height = 1.0;
    snapshot.total_mass = 2.0;

    Snapshot::Cell first;
    first.site_index = 0;
    first.site = Vector3(0.5, 0.5, 0.0);
    first.boundary = {Vector3(0.0, 0.0, 0.0), Vector3(1.0, 0.0, 0.0), Vector3(1.0, 1.0, 0.0), Vector3(0.0, 1.0, 0.0)};
    first.neighbors = {1};
    first.mass = 1.0;
    first.area = 1.0;

    Snapshot::Cell second = first;
    second.site_index = 1;
    second.site = Vector3(1.5, 0.5, 0.0);
    second.neighbors = {0};

    snapshot.cells = {first, second};

    return snapshot;
}

TEST(JsonReaderTest, ReadsWhatTheWriterWrote) {
    Snapshot written = two_cells();
    Snapshot read = JsonReader().read(JsonWriter().to_string(written));

    EXPECT_EQ(read.geometry, "plane");
    EXPECT_DOUBLE_EQ(read.width, 2.0);
    EXPECT_DOUBLE_EQ(read.height, 1.0);
    EXPECT_DOUBLE_EQ(read.total_mass, 2.0);
    ASSERT_EQ(read.cells.size(), 2u);

    for (size_t index = 0; index < 2; ++index) {
        EXPECT_EQ(read.cells[index].site_index, written.cells[index].site_index);
        EXPECT_NEAR((read.cells[index].site - written.cells[index].site).norm(), 0.0, 1e-9);
        EXPECT_EQ(read.cells[index].neighbors, written.cells[index].neighbors);
        EXPECT_EQ(read.cells[index].boundary.size(), 4u);
        EXPECT_DOUBLE_EQ(read.cells[index].mass, 1.0);
        EXPECT_DOUBLE_EQ(read.cells[index].area, 1.0);
    }

    EXPECT_NEAR((read.cells[1].boundary[2] - Vector3(1.0, 1.0, 0.0)).norm(), 0.0, 1e-9);
}

TEST(JsonReaderTest, ReadsAGoldenSphereSnapshot) {
    Snapshot snapshot = JsonReader().read_file(testdata("fibonacci_24.json"));

    EXPECT_EQ(snapshot.geometry, "sphere");
    ASSERT_EQ(snapshot.cells.size(), 24u);
    EXPECT_EQ(snapshot.cells[0].neighbors, (std::vector<size_t>{3, 8, 5, 2, 1}));
    EXPECT_EQ(snapshot.cells[0].boundary.size(), 5u);
    EXPECT_NEAR(snapshot.cells[0].site.norm(), 1.0, 1e-8);
}

TEST(JsonReaderTest, RefusesWhatIsNotASnapshot) {
    EXPECT_THROW((void) JsonReader().read("{\"particles\": []}"), std::runtime_error);
    EXPECT_THROW((void) JsonReader().read("not json"), std::runtime_error);
}

} // namespace
} // namespace geometry_art::io::snapshot
