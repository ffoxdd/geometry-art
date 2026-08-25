#ifndef GEOMETRY_ART_GENERATORS_SPHERICAL_POISSON_ELIMINATION_POINT_GENERATOR_HPP_
#define GEOMETRY_ART_GENERATORS_SPHERICAL_POISSON_ELIMINATION_POINT_GENERATOR_HPP_

#include "../../../geometry/spherical/indexed_kd_tree.hpp"
#include "../../../geometry/spherical/bounding_box.hpp"
#include "../../../geometry/spherical/helpers.hpp"
#include "../point_generator.hpp"
#include "../random_point_generator.hpp"
#include <boost/iterator/counting_iterator.hpp>
#include <vector>
#include <cmath>
#include <queue>
#include <algorithm>

namespace geometry_art::generators::spherical::poisson {

using geometry_art::geometry::spherical::IndexedPointMap;
using geometry_art::geometry::spherical::IndexedSearchTraits;
using geometry_art::geometry::spherical::IndexedKDTree;
using geometry_art::geometry::spherical::IndexedFuzzySphere;

template<spherical::PointGenerator PointGeneratorType = RandomPointGenerator<>>
class EliminationPointGenerator {
 public:
    EliminationPointGenerator() = default;

    explicit EliminationPointGenerator(
        PointGeneratorType generator,
        double oversample_factor = 2.0
    ) : _generator(std::move(generator)),
        _oversample_factor(oversample_factor) {
    }

    std::vector<VectorS2> generate(size_t count);
    std::vector<VectorS2> generate(size_t count, const SphericalBoundingBox &bounding_box);

    [[nodiscard]] size_t last_attempt_count() const { return _last_attempt_count; }

 private:
    PointGeneratorType _generator;
    double _oversample_factor = 2.0;
    size_t _last_attempt_count = 0;

    static constexpr double WEIGHT_ALPHA = 8.0;

    double compute_r_max(size_t target_count, double area) const;
    double weight_contribution(double geodesic_distance, double r_max) const;

    std::vector<VectorS2> eliminate_to_count(
        std::vector<VectorS2> candidates,
        size_t target_count,
        double r_max
    );
};

template<spherical::PointGenerator PointGeneratorType>
std::vector<VectorS2> EliminationPointGenerator<PointGeneratorType>::generate(size_t count) {
    return generate(count, SphericalBoundingBox::full_sphere());
}

template<spherical::PointGenerator PointGeneratorType>
std::vector<VectorS2> EliminationPointGenerator<PointGeneratorType>::generate(
    size_t count,
    const SphericalBoundingBox &bounding_box
) {
    if (count == 0) {
        _last_attempt_count = 0;
        return {};
    }

    size_t oversample_count = static_cast<size_t>(count * _oversample_factor);
    oversample_count = std::max(oversample_count, count + 1);

    auto candidates = _generator.generate(oversample_count, bounding_box);
    _last_attempt_count = _generator.last_attempt_count();

    double area = bounding_box.area();
    double r_max = compute_r_max(count, area);

    return eliminate_to_count(std::move(candidates), count, r_max);
}

template<spherical::PointGenerator PointGeneratorType>
double EliminationPointGenerator<PointGeneratorType>::compute_r_max(
    size_t target_count,
    double area
) const {
    // For hexagonal packing on a sphere, area per point ≈ (sqrt(3)/2) * r^2
    // Solving for r: r = sqrt(2 * area / (sqrt(3) * n))
    // We use 2x this as r_max for the neighborhood radius
    double area_per_point = area / static_cast<double>(target_count);
    double r_optimal = std::sqrt(2.0 * area_per_point / std::sqrt(3.0));
    return 2.0 * r_optimal;
}

template<spherical::PointGenerator PointGeneratorType>
double EliminationPointGenerator<PointGeneratorType>::weight_contribution(
    double geodesic_distance,
    double r_max
) const {
    if (geodesic_distance >= r_max) {
        return 0.0;
    }
    double ratio = 1.0 - geodesic_distance / r_max;
    return std::pow(ratio, WEIGHT_ALPHA);
}

template<spherical::PointGenerator PointGeneratorType>
std::vector<VectorS2> EliminationPointGenerator<PointGeneratorType>::eliminate_to_count(
    std::vector<VectorS2> candidates,
    size_t target_count,
    double r_max
) {
    size_t n = candidates.size();
    if (n <= target_count) {
        return candidates;
    }

    // Convert to cgal::Point3 for KD-tree (CGAL boundary)
    std::vector<cgal::Point3> cgal_points;
    cgal_points.reserve(n);
    for (const auto& p : candidates) {
        cgal_points.push_back(cgal::to_point(p));
    }

    // Chord distance corresponding to geodesic r_max: 2*sin(r_max/2)
    double chord_r_max = 2.0 * std::sin(r_max / 2.0);

    // Build indexed KD-tree - stores indices, maps to points via property map
    IndexedPointMap<std::vector<cgal::Point3>> point_map(cgal_points);
    IndexedSearchTraits search_traits(point_map);

    IndexedKDTree tree(
        boost::counting_iterator<std::size_t>(0),
        boost::counting_iterator<std::size_t>(n),
        IndexedKDTree::Splitter(),
        search_traits
    );
    tree.build();

    // Track which points are still active
    std::vector<bool> active(n, true);
    size_t active_count = n;

    // Precompute neighbor lists - queries now return indices directly
    std::vector<std::vector<size_t>> neighbor_indices(n);

    for (size_t i = 0; i < n; ++i) {
        IndexedFuzzySphere query(cgal_points[i], chord_r_max, 0.0, search_traits);
        std::vector<size_t> neighbors;
        tree.search(std::back_inserter(neighbors), query);

        for (size_t j : neighbors) {
            if (j != i) {
                neighbor_indices[i].push_back(j);
            }
        }
    }

    // Calculate initial weights using geodesic distance
    std::vector<double> weights(n, 0.0);

    for (size_t i = 0; i < n; ++i) {
        const VectorS2& vi = candidates[i];
        for (size_t j : neighbor_indices[i]) {
            const VectorS2& vj = candidates[j];
            double cos_theta = std::clamp(vi.dot(vj), -1.0, 1.0);
            double geodesic = std::acos(cos_theta);
            weights[i] += weight_contribution(geodesic, r_max);
        }
    }

    // Max-heap: (weight, index)
    using HeapEntry = std::pair<double, size_t>;
    std::priority_queue<HeapEntry> heap;

    for (size_t i = 0; i < n; ++i) {
        heap.emplace(weights[i], i);
    }

    // Eliminate points until we reach target count
    while (active_count > target_count && !heap.empty()) {
        auto [weight, index] = heap.top();
        heap.pop();

        if (!active[index]) {
            continue;
        }

        // Check if weight is stale (was updated after being pushed)
        if (std::abs(weight - weights[index]) > 1e-10) {
            // Re-push with current weight
            heap.emplace(weights[index], index);
            continue;
        }

        // Eliminate this point
        active[index] = false;
        active_count--;

        // Update weights of neighbors
        const VectorS2& vi = candidates[index];

        for (size_t j : neighbor_indices[index]) {
            if (!active[j]) continue;

            const VectorS2& vj = candidates[j];
            double cos_theta = std::clamp(vi.dot(vj), -1.0, 1.0);
            double geodesic = std::acos(cos_theta);
            double contribution = weight_contribution(geodesic, r_max);
            weights[j] -= contribution;
            // Push updated weight (lazy update - old entry will be ignored)
            heap.emplace(weights[j], j);
        }
    }

    // Collect remaining points
    std::vector<VectorS2> result;
    result.reserve(target_count);

    for (size_t i = 0; i < n && result.size() < target_count; ++i) {
        if (active[i]) {
            result.push_back(candidates[i]);
        }
    }

    return result;
}

} // namespace geometry_art::generators::spherical::poisson

#endif //GEOMETRY_ART_GENERATORS_SPHERICAL_POISSON_ELIMINATION_POINT_GENERATOR_HPP_
