// The top of the pyramid: the pipeline as a whole, and then the programs as
// programs. Nothing here reaches inside the solver -- these hold the
// properties a finished tessellation must have however it was reached.
#include "geometry_art/io/snapshot/snapshot.hpp"
#include "geometry_art/testing/macros.hpp"
#include "geometry_art/voronoi/spherical/factories/factory.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <set>
#include <sstream>
#include <algorithm>
#include <string>

using namespace geometry_art;
using io::snapshot::Snapshot;
using voronoi::CapacityConstrainedParameters;
using voronoi::spherical::Factory;
using voronoi::spherical::noop_callback;

namespace {

constexpr unsigned int SEED = 7;

std::filesystem::path scratch(const std::string& name) {
    std::filesystem::path directory = std::filesystem::temp_directory_path() / "geometry-art-pipeline-test";
    std::filesystem::create_directories(directory);
    return directory / name;
}

Snapshot built(int sites, const std::string& field, const std::string& inner_solver) {
    CapacityConstrainedParameters parameters;
    parameters.relative_capacity_tolerance = 1e-6;
    parameters.inner_solver = inner_solver;

    Factory factory(sites, field, 5, "lloyd", 0, parameters, SEED, noop_callback(), {}, std::chrono::milliseconds(0));
    factory.build();

    return factory.snapshot();
}

std::string contents_of(const std::filesystem::path& path) {
    std::ifstream stream(path);
    std::ostringstream buffer;
    buffer << stream.rdbuf();
    return buffer.str();
}

int run(const std::string& command) {
    return std::system((command + " > /dev/null 2>&1").c_str());
}

std::string binary(const std::string& name) {
    return std::string(GEOMETRY_ART_BINARY_DIR) + "/" + name;
}

} // namespace

TEST(PipelineTest, EveryCellIsAClosedPolygonWithNeighborsInRange) {
    Snapshot snapshot = built(40, "constant", "lbfgs");

    ASSERT_EQ(snapshot.cells.size(), 40u);

    for (const Snapshot::Cell& cell : snapshot.cells) {
        EXPECT_GE(cell.boundary.size(), 3u) << "cell " << cell.site_index;
        EXPECT_GT(cell.area, 0.0) << "cell " << cell.site_index;
        EXPECT_NEAR(cell.site.norm(), 1.0, 1e-9);

        for (size_t neighbor : cell.neighbors) {
            EXPECT_LT(neighbor, snapshot.cells.size());
            EXPECT_NE(neighbor, cell.site_index);
        }

        for (const Vector3& point : cell.boundary) {
            EXPECT_NEAR(point.norm(), 1.0, 1e-9);
        }
    }
}

TEST(PipelineTest, CellsPartitionTheSphere) {
    Snapshot snapshot = built(40, "constant", "lbfgs");
    double area = 0.0;
    double mass = 0.0;

    for (const Snapshot::Cell& cell : snapshot.cells) {
        area += cell.area;
        mass += cell.mass;
    }

    EXPECT_NEAR(area, 4.0 * M_PI, 1e-9);
    EXPECT_NEAR(mass, snapshot.total_mass, 1e-9);
}

TEST(PipelineTest, AdjacencyIsSymmetric) {
    Snapshot snapshot = built(40, "constant", "lbfgs");

    for (const Snapshot::Cell& cell : snapshot.cells) {
        for (size_t neighbor : cell.neighbors) {
            const std::vector<size_t>& theirs = snapshot.cells[neighbor].neighbors;
            EXPECT_NE(std::find(theirs.begin(), theirs.end(), cell.site_index), theirs.end())
                << cell.site_index << " lists " << neighbor << " but not the reverse";
        }
    }
}

TEST(PipelineTest, CapacitiesAreEqualisedToTheRequestedTolerance) {
    Snapshot snapshot = built(40, "quadratic", "lbfgs");

    EXPECT_LT(snapshot.relative_rms_capacity_error(), 1e-6);
}

TEST(PipelineTest, BothInnerSolversReachTheSameTessellationQuality) {
    Snapshot lbfgs = built(40, "quadratic", "lbfgs");
    Snapshot newton = built(40, "quadratic", "newton");

    EXPECT_LT(newton.relative_rms_capacity_error(), 1e-6);
    EXPECT_NEAR(lbfgs.total_mass, newton.total_mass, 1e-12);
    EXPECT_EQ(lbfgs.cells.size(), newton.cells.size());
}

TEST(PipelineTest, EXPENSIVE_TessellateWritesAReadableSnapshot) {
    REQUIRE_EXPENSIVE();

    std::filesystem::path prefix = scratch("cli");
    std::filesystem::remove(prefix.string() + ".json");
    std::filesystem::remove(prefix.string() + ".svg");

    int status = run(
        binary("tessellate") + " -f constant -p 24 --seed 7"
        " --capacity-tolerance 1e-6"
        " --snapshot " + prefix.string() + " -o " + scratch("out").string()
    );

    ASSERT_EQ(status, 0);
    ASSERT_TRUE(std::filesystem::exists(prefix.string() + ".json"));
    ASSERT_TRUE(std::filesystem::exists(prefix.string() + ".svg"));

    std::string json = contents_of(prefix.string() + ".json");
    EXPECT_NE(json.find("\"geometry\": \"sphere\""), std::string::npos);
    EXPECT_NE(json.find("\"relativeRmsCapacityError\""), std::string::npos);

    size_t cells = 0;
    size_t position = json.find("{\"site\":");

    while (position != std::string::npos) {
        ++cells;
        position = json.find("{\"site\":", position + 1);
    }

    EXPECT_EQ(cells, 24u);
    EXPECT_NE(contents_of(prefix.string() + ".svg").find("<svg"), std::string::npos);
}

TEST(PipelineTest, EXPENSIVE_SkeletonizeWritesAClosedSolidFromASnapshot) {
    REQUIRE_EXPENSIVE();

    std::filesystem::path prefix = scratch("skeleton");
    std::filesystem::remove(prefix.string() + ".json");
    std::filesystem::remove(prefix.string() + ".stl");
    std::filesystem::remove(prefix.string() + ".obj");

    ASSERT_EQ(run(
        binary("tessellate") + " -f constant -p 24 --seed 7"
        " --capacity-tolerance 1e-6 --snapshot " + prefix.string()
    ), 0);

    // The model is named after the snapshot unless told otherwise, and the
    // default format is STL, written in its binary form: an 80 byte header,
    // a facet count, and fifty bytes per facet.
    ASSERT_EQ(run(binary("skeletonize") + " " + prefix.string() + ".json --resolution 2.0"), 0);

    std::filesystem::path stl = prefix.string() + ".stl";
    ASSERT_TRUE(std::filesystem::exists(stl));

    std::string body = contents_of(stl);
    ASSERT_GE(body.size(), 84u);

    uint32_t facets = 0;
    std::memcpy(&facets, body.data() + 80, sizeof(facets));

    EXPECT_GT(facets, 0u);
    EXPECT_EQ(body.size(), 84u + 50u * facets);

    ASSERT_EQ(run(binary("skeletonize") + " " + prefix.string() + ".json -f obj --resolution 2.0"), 0);

    std::string obj = contents_of(prefix.string() + ".obj");
    EXPECT_NE(obj.find("\nv "), std::string::npos);
    EXPECT_NE(obj.find("\nf "), std::string::npos);
}
