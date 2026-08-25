#ifndef GEOMETRY_ART_VORONOI_SPHERICAL_OPTIMIZERS_LLOYD_OPTIMIZER_HPP_
#define GEOMETRY_ART_VORONOI_SPHERICAL_OPTIMIZERS_LLOYD_OPTIMIZER_HPP_

#include "../../../types.hpp"
#include "../../../fields/spherical/field.hpp"
#include "../../../fields/spherical/polynomial_field.hpp"
#include "../../../geometry/spherical/helpers.hpp"
#include "../core/sphere.hpp"
#include "../core/callback.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <memory>
#include <utility>

namespace geometry_art::voronoi::spherical {

using geometry::spherical::distance;

template<fields::spherical::Field FieldType = fields::spherical::PolynomialField>
class LloydOptimizer {
 public:
    LloydOptimizer(
        std::unique_ptr<Sphere> sphere,
        FieldType field,
        size_t passes,
        Callback callback
    );

    std::unique_ptr<Sphere> optimize();
    [[nodiscard]] double final_deviation() const { return _final_deviation; }

 private:
    std::unique_ptr<Sphere> _sphere;
    FieldType _field;
    size_t _passes;
    Callback _callback;
    double _final_deviation = 0.0;

    void run_single_pass();
    [[nodiscard]] double root_mean_square_deviation() const;
    [[nodiscard]] VectorS2 weighted_centroid(const Polygon& cell) const;
};

template<fields::spherical::Field FieldType>
LloydOptimizer<FieldType>::LloydOptimizer(
    std::unique_ptr<Sphere> sphere,
    FieldType field,
    size_t passes,
    Callback callback
) :
    _sphere(std::move(sphere)),
    _field(std::move(field)),
    _passes(passes),
    _callback(std::move(callback)) {
}

template<fields::spherical::Field FieldType>
std::unique_ptr<Sphere> LloydOptimizer<FieldType>::optimize() {
    for (size_t pass = 0; pass < _passes; ++pass) {
        run_single_pass();
        _callback(*_sphere);
    }

    _final_deviation = root_mean_square_deviation();
    return std::move(_sphere);
}

template<fields::spherical::Field FieldType>
void LloydOptimizer<FieldType>::run_single_pass() {
    std::vector<VectorS2> centroids;
    centroids.reserve(_sphere->size());

    for (const auto& cell : _sphere->cells()) {
        centroids.push_back(weighted_centroid(cell));
    }

    for (size_t index = 0; index < centroids.size(); ++index) {
        _sphere->update_site(index, cgal::to_point(centroids[index]));
    }
}

template<fields::spherical::Field FieldType>
double LloydOptimizer<FieldType>::root_mean_square_deviation() const {
    double total = 0.0;
    size_t index = 0;

    for (const auto& cell : _sphere->cells()) {
        double deviation = distance(to_vector_s2(_sphere->site(index)), weighted_centroid(cell));
        total += deviation * deviation;
        ++index;
    }

    return std::sqrt(total / static_cast<double>(_sphere->size()));
}

template<fields::spherical::Field FieldType>
VectorS2 LloydOptimizer<FieldType>::weighted_centroid(const Polygon& cell) const {
    Vector3 moment = _field.integrals(cell).first_moment;
    double norm = moment.norm();
    return norm < GEOMETRIC_EPSILON ? cell.centroid() : VectorS2(moment / norm);
}

} // namespace geometry_art::voronoi::spherical

#endif //GEOMETRY_ART_VORONOI_SPHERICAL_OPTIMIZERS_LLOYD_OPTIMIZER_HPP_
