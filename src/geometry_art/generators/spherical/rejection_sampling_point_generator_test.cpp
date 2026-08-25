#include <gtest/gtest.h>
#include "rejection_sampling_point_generator.hpp"
#include "../../fields/scalar/constant_field.hpp"
#include "../../testing/assertions/geometric.hpp"

using namespace geometry_art::generators::spherical;
using geometry_art::fields::scalar::ConstantField;
using geometry_art::VectorS2;
using geometry_art::testing::is_on_unit_sphere;

TEST(RejectionSamplingPointGeneratorTest, GenerateReturnsExactCount) {
    RejectionSamplingPointGenerator generator;
    constexpr size_t count = 100;

    auto points = generator.generate(count);

    EXPECT_EQ(points.size(), count);
}

TEST(RejectionSamplingPointGeneratorTest, GenerateZeroReturnsEmpty) {
    RejectionSamplingPointGenerator generator;

    auto points = generator.generate(0);

    EXPECT_TRUE(points.empty());
    EXPECT_EQ(generator.last_attempt_count(), 0);
}

TEST(RejectionSamplingPointGeneratorTest, AllPointsOnUnitSphere) {
    RejectionSamplingPointGenerator generator;
    constexpr size_t count = 50;

    auto points = generator.generate(count);

    for (const auto &point : points) {
        EXPECT_TRUE(is_on_unit_sphere(point));
    }
}

TEST(RejectionSamplingPointGeneratorTest, UniformDensityAcceptsAll) {
    ConstantField uniform_field(1.0);
    RejectionSamplingPointGenerator generator(uniform_field, 1.0);
    constexpr size_t count = 100;

    auto points = generator.generate(count);

    EXPECT_EQ(points.size(), count);
    EXPECT_EQ(generator.last_attempt_count(), count);
}

TEST(RejectionSamplingPointGeneratorTest, HalfDensityRejectsApproximatelyHalf) {
    ConstantField half_density_field(0.5);
    RejectionSamplingPointGenerator generator(half_density_field, 1.0);
    constexpr size_t count = 100;

    auto points = generator.generate(count);

    EXPECT_EQ(points.size(), count);
    EXPECT_GT(generator.last_attempt_count(), count);
    EXPECT_LT(generator.last_attempt_count(), count * 4);
}

TEST(RejectionSamplingPointGeneratorTest, ZeroDensityNeverAccepts) {
    ConstantField zero_field(0.0);
    RejectionSamplingPointGenerator generator(zero_field, 1.0);

    auto points = generator.generate(0);

    EXPECT_TRUE(points.empty());
}

TEST(RejectionSamplingPointGeneratorTest, LastAttemptCountTracksAttempts) {
    ConstantField field(0.25);
    RejectionSamplingPointGenerator generator(field, 1.0);
    constexpr size_t count = 50;

    generator.generate(count);

    EXPECT_GE(generator.last_attempt_count(), count);
}
