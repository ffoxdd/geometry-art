#include <gtest/gtest.h>
#include "constant_field.hpp"

using namespace geometry_art::fields::scalar;
using geometry_art::VectorS2;

TEST(ConstantFieldTest, ReturnsConfiguredValue) {
    ConstantField field(1.0);
    EXPECT_DOUBLE_EQ(field.value(VectorS2(0, 0, 1)), 1.0);

    ConstantField other_field(2.5);
    EXPECT_DOUBLE_EQ(other_field.value(VectorS2(1, 0, 0)), 2.5);
}

