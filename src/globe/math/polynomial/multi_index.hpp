#ifndef GLOBEART_SRC_GLOBE_MATH_POLYNOMIAL_MULTI_INDEX_HPP_
#define GLOBEART_SRC_GLOBE_MATH_POLYNOMIAL_MULTI_INDEX_HPP_

#include "maximum_degree.hpp"
#include <CGAL/assertions.h>
#include <array>
#include <cstddef>
#include <vector>

namespace globe::math::polynomial {

struct MultiIndex {
    int x = 0;
    int y = 0;
    int z = 0;

    [[nodiscard]] constexpr int degree() const { return x + y + z; }
    [[nodiscard]] constexpr int operator[](int axis) const { return axis == 0 ? x : (axis == 1 ? y : z); }
    [[nodiscard]] constexpr MultiIndex raised(int axis) const;
    [[nodiscard]] constexpr MultiIndex lowered(int axis) const;
    [[nodiscard]] constexpr bool operator==(const MultiIndex& other) const = default;

    [[nodiscard]] static constexpr MultiIndex unit(int axis);
    [[nodiscard]] static constexpr size_t count_up_to(int max_degree);
    [[nodiscard]] static constexpr size_t ordinal(const MultiIndex& index);
    [[nodiscard]] static const std::vector<MultiIndex>& all_up_to(int max_degree);
    [[nodiscard]] static const std::vector<MultiIndex>& all_of_degree(int degree);
};

constexpr MultiIndex MultiIndex::unit(int axis) {
    return MultiIndex{axis == 0 ? 1 : 0, axis == 1 ? 1 : 0, axis == 2 ? 1 : 0};
}

constexpr MultiIndex MultiIndex::raised(int axis) const {
    return MultiIndex{x + (axis == 0), y + (axis == 1), z + (axis == 2)};
}

constexpr MultiIndex MultiIndex::lowered(int axis) const {
    return MultiIndex{x - (axis == 0), y - (axis == 1), z - (axis == 2)};
}

constexpr size_t MultiIndex::count_up_to(int max_degree) {
    size_t d = static_cast<size_t>(max_degree);
    return (d + 1) * (d + 2) * (d + 3) / 6;
}

constexpr size_t MultiIndex::ordinal(const MultiIndex& index) {
    size_t d = static_cast<size_t>(index.degree());
    size_t degree_block_offset = d * (d + 1) * (d + 2) / 6;
    size_t remainder = static_cast<size_t>(index.y + index.z);
    return degree_block_offset + remainder * (remainder + 1) / 2 + static_cast<size_t>(index.z);
}

// The enumerations depend on nothing but the degree, and the moment
// recursions walk them once per region, so each is built once and shared.
inline const std::vector<MultiIndex>& MultiIndex::all_up_to(int max_degree) {
    CGAL_precondition(max_degree >= 0 && max_degree <= MAXIMUM_DEGREE);

    static const std::vector<std::vector<MultiIndex>> table = [] {
        std::vector<std::vector<MultiIndex>> result(MAXIMUM_DEGREE + 1);

        for (int max = 0; max <= MAXIMUM_DEGREE; ++max) {
            result[max].reserve(count_up_to(max));

            for (int degree = 0; degree <= max; ++degree) {
                for (const MultiIndex& index : all_of_degree(degree)) {
                    result[max].push_back(index);
                }
            }
        }

        return result;
    }();

    return table[max_degree];
}

inline const std::vector<MultiIndex>& MultiIndex::all_of_degree(int degree) {
    CGAL_precondition(degree >= 0 && degree <= MAXIMUM_DEGREE);

    static const std::vector<std::vector<MultiIndex>> table = [] {
        std::vector<std::vector<MultiIndex>> result(MAXIMUM_DEGREE + 1);

        for (int total = 0; total <= MAXIMUM_DEGREE; ++total) {
            for (int x = total; x >= 0; --x) {
                for (int y = total - x; y >= 0; --y) {
                    result[total].push_back(MultiIndex{x, y, total - x - y});
                }
            }
        }

        return result;
    }();

    return table[degree];
}

} // namespace globe::math::polynomial

#endif //GLOBEART_SRC_GLOBE_MATH_POLYNOMIAL_MULTI_INDEX_HPP_
