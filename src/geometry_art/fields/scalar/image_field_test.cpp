#include "image_field.hpp"
#include "../spherical/powell_sabin_projection.hpp"
#include "../../geometry/spherical/triangle_mesh.hpp"
#include "../../generators/spherical/fibonacci_point_generator.hpp"
#include "../../math/interval.hpp"
#include "../../testing/macros.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <vector>

using namespace geometry_art;
using fields::scalar::ImageField;
using fields::spherical::PowellSabinProjection;
using geometry::spherical::TriangleMesh;

namespace {

ImageField uniform_image(double brightness, Interval range) {
    return ImageField(std::vector<double>(8, brightness), 4, 2, range);
}

// Black on the western hemisphere, white on the eastern: a step in
// longitude, the hardest input a continuous density can be asked to hold.
ImageField step_image(size_t width, size_t height, Interval range) {
    std::vector<double> brightness;
    brightness.reserve(width * height);

    for (size_t row = 0; row < height; ++row) {
        for (size_t column = 0; column < width; ++column) {
            brightness.push_back(column < width / 2 ? 0.0 : 1.0);
        }
    }

    return ImageField(std::move(brightness), width, height, range);
}

} // namespace

// Darker is denser: black reads as the top of the range, white as the
// floor, so the cells crowd where the ink is.
TEST(ImageFieldTest, BlackIsDenseAndWhiteIsTheFloor) {
    Interval range(0.2, 1.0);
    VectorS2 point = VectorS2(1.0, 0.3, -0.2).normalized();

    EXPECT_NEAR(uniform_image(0.0, range).value(point), 1.0, 1e-12);
    EXPECT_NEAR(uniform_image(1.0, range).value(point), 0.2, 1e-12);
    EXPECT_NEAR(uniform_image(0.5, range).value(point), 0.6, 1e-12);
}

TEST(ImageFieldTest, IsContinuousAcrossTheLongitudeSeam) {
    ImageField image = step_image(16, 8, Interval(0.2, 1.0));

    double just_west = image.value(VectorS2(std::cos(M_PI - 1e-6), std::sin(M_PI - 1e-6), 0.0));
    double just_east = image.value(VectorS2(std::cos(M_PI + 1e-6), std::sin(M_PI + 1e-6), 0.0));

    EXPECT_NEAR(just_west, just_east, 1e-4);
}

TEST(ImageFieldTest, InterpolatesBetweenPixelCenters) {
    // One black and one white pixel around the equator: halfway between
    // their centers the density is the middle of the range.
    ImageField image(std::vector<double>{0.0, 1.0}, 2, 1, Interval(0.0, 1.0));

    EXPECT_NEAR(image.value(VectorS2(1.0, 0.0, 0.0)), 0.5, 1e-12);
    EXPECT_NEAR(image.value(VectorS2(-1.0, 0.0, 0.0)), 0.5, 1e-6);
}

// The verification the image path is held to: a step becomes a ramp about
// one cell wide, and the projection does not ring below the floor on its
// way down. A C1 field cannot hold the step, and the bandwidth rule says
// the tessellation could not read it anyway.
TEST(ImageFieldTest, EXPENSIVE_ProjectsAStepEdgeWithoutRinging) {
    REQUIRE_EXPENSIVE();

    Interval range(0.2, 1.0);
    ImageField image = step_image(64, 32, range);
    PowellSabinProjection projection;
    auto result = projection.project(TriangleMesh::icosphere(3), image);

    EXPECT_GT(result.lowest_coefficient, 0.0);

    double lowest = 1e9;
    double highest = -1e9;

    for (const VectorS2& point : generators::spherical::FibonacciPointGenerator().generate(20000)) {
        double value = result.field.value(point);
        lowest = std::min(lowest, value);
        highest = std::max(highest, value);
    }

    EXPECT_GT(lowest, range.low() - 0.1 * range.measure());
    EXPECT_LT(highest, range.high() + 0.1 * range.measure());

    // Far from the edge the plateaus are held exactly.
    EXPECT_NEAR(result.field.value(VectorS2(0.0, -1.0, 0.0).normalized()), range.high(), 0.02);
    EXPECT_NEAR(result.field.value(VectorS2(0.0, 1.0, 0.0).normalized()), range.low(), 0.02);
}
