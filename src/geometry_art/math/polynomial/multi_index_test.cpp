#include "multi_index.hpp"
#include <gtest/gtest.h>
#include <set>

using namespace geometry_art::math::polynomial;

TEST(MultiIndexTest, DegreeIsSumOfExponents) {
    EXPECT_EQ((MultiIndex{1, 2, 3}).degree(), 6);
}

TEST(MultiIndexTest, RaisedAndLoweredAdjustOneAxis) {
    MultiIndex index{1, 1, 1};

    EXPECT_EQ(index.raised(0), (MultiIndex{2, 1, 1}));
    EXPECT_EQ(index.raised(2), (MultiIndex{1, 1, 2}));
    EXPECT_EQ(index.lowered(1), (MultiIndex{1, 0, 1}));
}

TEST(MultiIndexTest, CountUpToMatchesTetrahedralNumbers) {
    EXPECT_EQ(MultiIndex::count_up_to(0), 1u);
    EXPECT_EQ(MultiIndex::count_up_to(1), 4u);
    EXPECT_EQ(MultiIndex::count_up_to(2), 10u);
    EXPECT_EQ(MultiIndex::count_up_to(3), 20u);
}

TEST(MultiIndexTest, OrdinalIsABijectionOntoDenseRange) {
    constexpr int MAX_DEGREE = 5;
    std::set<size_t> ordinals;

    for (const MultiIndex& index : MultiIndex::all_up_to(MAX_DEGREE)) {
        ordinals.insert(MultiIndex::ordinal(index));
    }

    EXPECT_EQ(ordinals.size(), MultiIndex::count_up_to(MAX_DEGREE));
    EXPECT_EQ(*ordinals.begin(), 0u);
    EXPECT_EQ(*ordinals.rbegin(), MultiIndex::count_up_to(MAX_DEGREE) - 1);
}

TEST(MultiIndexTest, OrdinalFollowsEnumerationOrder) {
    size_t expected = 0;

    for (const MultiIndex& index : MultiIndex::all_up_to(4)) {
        EXPECT_EQ(MultiIndex::ordinal(index), expected);
        ++expected;
    }
}

TEST(MultiIndexTest, AllOfDegreeListsEveryCombination) {
    auto indices = MultiIndex::all_of_degree(2);

    ASSERT_EQ(indices.size(), 6u);
    EXPECT_EQ(indices.front(), (MultiIndex{2, 0, 0}));
    EXPECT_EQ(indices.back(), (MultiIndex{0, 0, 2}));
}
