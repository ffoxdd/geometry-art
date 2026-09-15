#ifndef GEOMETRY_ART_SKELETON_ANNULUS_HPP_
#define GEOMETRY_ART_SKELETON_ANNULUS_HPP_

#include "../types.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

namespace geometry_art::skeleton {

using Triangle = std::array<size_t, 3>;

// Triangulates the ring between two nested loops, both counter-clockwise
// and star-shaped about the inner loop's centroid, by walking the two loops
// together in angle: whichever loop's next vertex comes first closes a
// triangle with the current vertex of the other. Triangles come out
// counter-clockwise, and their indices count the outer loop first. Without
// an inner loop the outer one is fanned from its first vertex.
[[nodiscard]] std::vector<Triangle> triangulate_annulus(
    const std::vector<Vector2>& outer,
    const std::vector<Vector2>& inner
);

namespace annulus {

// The vertices of a loop in travel order from its first vertex past the
// centre's leftmost ray, each with its angle about the centre unwrapped so
// the sequence climbs monotonically, and the first vertex repeated a full
// turn later to close the walk.
struct Walk {
    std::vector<size_t> indices;
    std::vector<double> angles;
};

[[nodiscard]] inline Walk walk_of(const std::vector<Vector2>& loop, const Vector2& center) {
    std::vector<double> angles;
    angles.reserve(loop.size());

    for (const Vector2& point : loop) {
        angles.push_back(std::atan2(point.y() - center.y(), point.x() - center.x()));
    }

    size_t start = static_cast<size_t>(std::min_element(angles.begin(), angles.end()) - angles.begin());
    Walk walk;
    double lift = 0.0;
    double previous = angles[start];

    for (size_t step = 0; step <= loop.size(); ++step) {
        size_t index = (start + step) % loop.size();
        double angle = angles[index] + (step == loop.size() ? TWO_PI : 0.0);

        if (angle < previous) {
            lift += TWO_PI;
        }

        walk.indices.push_back(index);
        walk.angles.push_back(angle + lift);
        previous = angle;
    }

    return walk;
}

[[nodiscard]] inline Vector2 centroid_of(const std::vector<Vector2>& loop) {
    Vector2 sum = Vector2::Zero();

    for (const Vector2& point : loop) {
        sum += point;
    }

    return sum / static_cast<double>(loop.size());
}

[[nodiscard]] inline std::vector<Triangle> fan(size_t count) {
    std::vector<Triangle> triangles;

    for (size_t index = 1; index + 1 < count; ++index) {
        triangles.push_back(Triangle{0, index, index + 1});
    }

    return triangles;
}

} // namespace annulus

inline std::vector<Triangle> triangulate_annulus(
    const std::vector<Vector2>& outer,
    const std::vector<Vector2>& inner
) {
    if (outer.size() < 3) {
        return {};
    }

    if (inner.size() < 3) {
        return annulus::fan(outer.size());
    }

    Vector2 center = annulus::centroid_of(inner);
    annulus::Walk outer_walk = annulus::walk_of(outer, center);
    annulus::Walk inner_walk = annulus::walk_of(inner, center);
    size_t outer_count = outer.size();
    size_t inner_count = inner.size();
    std::vector<Triangle> triangles;
    triangles.reserve(outer_count + inner_count);

    size_t outer_at = 0;
    size_t inner_at = 0;

    while (outer_at < outer_count || inner_at < inner_count) {
        bool advance_outer =
            inner_at == inner_count ||
            (outer_at < outer_count && outer_walk.angles[outer_at + 1] <= inner_walk.angles[inner_at + 1]);

        size_t outer_here = outer_walk.indices[outer_at];
        size_t inner_here = outer_count + inner_walk.indices[inner_at];

        if (advance_outer) {
            triangles.push_back(Triangle{outer_here, outer_walk.indices[outer_at + 1], inner_here});
            ++outer_at;
        } else {
            triangles.push_back(Triangle{outer_here, outer_count + inner_walk.indices[inner_at + 1], inner_here});
            ++inner_at;
        }
    }

    return triangles;
}

} // namespace geometry_art::skeleton

#endif //GEOMETRY_ART_SKELETON_ANNULUS_HPP_
