#ifndef GLOBEART_SRC_GLOBE_MATH_POLYNOMIAL_MULTI_INDEX_HPP_
#define GLOBEART_SRC_GLOBE_MATH_POLYNOMIAL_MULTI_INDEX_HPP_

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
    [[nodiscard]] static std::vector<MultiIndex> all_up_to(int max_degree);
    [[nodiscard]] static std::vector<MultiIndex> all_of_degree(int degree);
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

inline std::vector<MultiIndex> MultiIndex::all_up_to(int max_degree) {
    std::vector<MultiIndex> result;
    result.reserve(count_up_to(max_degree));

    for (int degree = 0; degree <= max_degree; ++degree) {
        for (const MultiIndex& index : all_of_degree(degree)) {
            result.push_back(index);
        }
    }

    return result;
}

inline std::vector<MultiIndex> MultiIndex::all_of_degree(int degree) {
    std::vector<MultiIndex> result;

    for (int x = degree; x >= 0; --x) {
        for (int y = degree - x; y >= 0; --y) {
            result.push_back(MultiIndex{x, y, degree - x - y});
        }
    }

    return result;
}

} // namespace globe::math::polynomial

#endif //GLOBEART_SRC_GLOBE_MATH_POLYNOMIAL_MULTI_INDEX_HPP_
