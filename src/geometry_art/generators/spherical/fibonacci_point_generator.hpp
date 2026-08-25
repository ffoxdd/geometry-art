#ifndef GEOMETRY_ART_GENERATORS_SPHERICAL_FIBONACCI_POINT_GENERATOR_HPP_
#define GEOMETRY_ART_GENERATORS_SPHERICAL_FIBONACCI_POINT_GENERATOR_HPP_

#include "point_generator.hpp"
#include "../../types.hpp"
#include "../../geometry/spherical/bounding_box.hpp"
#include <vector>
#include <cmath>

namespace geometry_art::generators::spherical {

class FibonacciPointGenerator {
 public:
    FibonacciPointGenerator() = default;

    std::vector<VectorS2> generate(size_t count);
    std::vector<VectorS2> generate(size_t count, const SphericalBoundingBox& bounding_box);

    [[nodiscard]] size_t last_attempt_count() const { return _last_attempt_count; }

 private:
    static constexpr double GOLDEN_RATIO = 1.618033988749895;

    size_t _last_attempt_count = 0;
};

inline std::vector<VectorS2> FibonacciPointGenerator::generate(size_t count) {
    std::vector<VectorS2> points;
    points.reserve(count);

    for (size_t i = 0; i < count; ++i) {
        double theta = 2.0 * M_PI * i / GOLDEN_RATIO;
        double z = 1.0 - 2.0 * (i + 0.5) / count;
        double r = std::sqrt(1.0 - z * z);

        points.emplace_back(r * std::cos(theta), r * std::sin(theta), z);
    }

    _last_attempt_count = count;
    return points;
}

inline std::vector<VectorS2> FibonacciPointGenerator::generate(
    size_t count,
    const SphericalBoundingBox& bounding_box
) {
    auto all_points = generate(count);

    std::vector<VectorS2> filtered;
    filtered.reserve(count);

    for (const auto& point : all_points) {
        if (bounding_box.contains(point)) {
            filtered.push_back(point);
        }
    }

    _last_attempt_count = count;
    return filtered;
}

static_assert(PointGenerator<FibonacciPointGenerator>);

} // namespace geometry_art::generators::spherical

#endif //GEOMETRY_ART_GENERATORS_SPHERICAL_FIBONACCI_POINT_GENERATOR_HPP_
