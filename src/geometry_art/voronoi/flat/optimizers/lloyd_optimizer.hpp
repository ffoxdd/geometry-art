#ifndef GEOMETRY_ART_VORONOI_FLAT_OPTIMIZERS_LLOYD_OPTIMIZER_HPP_
#define GEOMETRY_ART_VORONOI_FLAT_OPTIMIZERS_LLOYD_OPTIMIZER_HPP_

#include "../core/diagram.hpp"
#include "../../../fields/flat/field.hpp"
#include "../../../types.hpp"
#include <cmath>
#include <cstddef>
#include <functional>
#include <memory>
#include <utility>
#include <vector>

namespace geometry_art::voronoi::flat {

using Callback = std::function<void(const Diagram&)>;

inline Callback noop_callback() {
    return [](const Diagram&) {};
}

template<fields::flat::Field FieldType>
class LloydOptimizer {
 public:
    LloydOptimizer(std::unique_ptr<Diagram> diagram, FieldType field, size_t passes, Callback callback);

    std::unique_ptr<Diagram> optimize();
    [[nodiscard]] double final_deviation() const { return _final_deviation; }

 private:
    std::unique_ptr<Diagram> _torus;
    FieldType _field;
    size_t _passes;
    Callback _callback;
    double _final_deviation = 0.0;

    [[nodiscard]] std::vector<Vector3> weighted_centroids() const;
};

template<fields::flat::Field FieldType>
LloydOptimizer<FieldType>::LloydOptimizer(
    std::unique_ptr<Diagram> diagram,
    FieldType field,
    size_t passes,
    Callback callback
) :
    _torus(std::move(diagram)),
    _field(std::move(field)),
    _passes(passes),
    _callback(std::move(callback)) {
}

template<fields::flat::Field FieldType>
std::unique_ptr<Diagram> LloydOptimizer<FieldType>::optimize() {
    for (size_t pass = 0; pass < _passes; ++pass) {
        _torus = _torus->rebuilt(weighted_centroids());
        _callback(*_torus);
    }

    std::vector<Vector3> centroids = weighted_centroids();
    double total = 0.0;

    for (size_t index = 0; index < _torus->size(); ++index) {
        total += (centroids[index] - _torus->site_vector(index)).squaredNorm();
    }

    _final_deviation = std::sqrt(total / static_cast<double>(_torus->size()));
    return std::move(_torus);
}

// The centroid is taken in the cell's own chart, where the cell is a plain
// polygon; the rebuild wraps it back into the rectangle.
template<fields::flat::Field FieldType>
std::vector<Vector3> LloydOptimizer<FieldType>::weighted_centroids() const {
    std::vector<Vector3> centroids;
    centroids.reserve(_torus->size());

    for (size_t index = 0; index < _torus->size(); ++index) {
        auto integrals = _field.integrals(_torus->cell(index));

        centroids.push_back(
            integrals.mass < GEOMETRIC_EPSILON
                ? _torus->site_vector(index)
                : Vector3(integrals.first_moment / integrals.mass)
        );
    }

    return centroids;
}

} // namespace geometry_art::voronoi::flat

#endif //GEOMETRY_ART_VORONOI_FLAT_OPTIMIZERS_LLOYD_OPTIMIZER_HPP_
