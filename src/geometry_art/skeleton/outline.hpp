#ifndef GEOMETRY_ART_SKELETON_OUTLINE_HPP_
#define GEOMETRY_ART_SKELETON_OUTLINE_HPP_

#include "../types.hpp"
#include <cmath>
#include <concepts>
#include <cstddef>
#include <vector>

namespace geometry_art::skeleton {

// One point of a sampled boundary loop: where it sits in the cell's chart,
// which is what the loop is triangulated in, and where it sits in space
// together with the surface normal there. On an outer loop a sample also
// says whether the edge leaving it lies on the domain's boundary, since
// that edge carries a wall.
struct Sample {
    Vector2 chart;
    Vector3 position;
    Vector3 normal;
    bool wall = false;
};

using Loop = std::vector<Sample>;

// A cell's share of the skeleton: the ring between its boundary and its
// inset, both counter-clockwise in the chart. An empty inner loop is a cell
// the bars fill entirely.
struct Outline {
    Loop outer;
    Loop inner;
};

// A geometry that can outline each of its cells with bars of a given half
// width, sampling curved edges no coarser than a step.
template<typename T>
concept Outliner = requires(const T& outliner, size_t index, double half_width, double step) {
    { outliner.size() } -> std::convertible_to<size_t>;
    { outliner.outline(index, half_width, step) } -> std::same_as<Outline>;
};

// How many pieces an edge of a given length is sampled in.
[[nodiscard]] inline size_t segments(double length, double step) {
    if (!(step > 0.0) || !(length > 0.0)) {
        return 1;
    }

    return static_cast<size_t>(std::max(1.0, std::ceil(length / step)));
}

} // namespace geometry_art::skeleton

#endif //GEOMETRY_ART_SKELETON_OUTLINE_HPP_
