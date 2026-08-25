#ifndef GLOBEART_SRC_GLOBE_GEOMETRY_PLANAR_SEGMENT_HPP_
#define GLOBEART_SRC_GLOBE_GEOMETRY_PLANAR_SEGMENT_HPP_

#include "../../types.hpp"
#include "../../math/binomials.hpp"
#include "../../math/powers.hpp"
#include "../../math/polynomial/moments.hpp"
#include "../../math/polynomial/multi_index.hpp"
#include <cmath>
#include <vector>

namespace globe::geometry::planar {

using globe::math::polynomial::Moments;
using globe::math::polynomial::MultiIndex;
using globe::math::powers;

// A straight edge in the plane. Moments are reported in the same
// three-index table the spherical types use, with the third exponent left at
// zero, so the polynomial and field machinery is shared between geometries.
class Segment {
 public:
    Segment(const Vector2& source, const Vector2& target);

    [[nodiscard]] const Vector2& source() const { return _source; }
    [[nodiscard]] const Vector2& target() const { return _target; }
    [[nodiscard]] Vector2 direction() const { return _target - _source; }
    [[nodiscard]] double length() const { return direction().norm(); }
    [[nodiscard]] Vector2 interpolate(double fraction) const { return _source + fraction * direction(); }

    // The outward normal for a counter-clockwise boundary.
    [[nodiscard]] Vector2 normal() const;

    [[nodiscard]] Moments moments(int max_degree) const;

    // Integrals of each monomial along the parametrisation, without the
    // length factor: the polygon's own moments need these weighted by a
    // component of the direction rather than by arc length.
    [[nodiscard]] Moments parametric_moments(int max_degree) const;

 private:
    Vector2 _source;
    Vector2 _target;
};

inline Segment::Segment(const Vector2& source, const Vector2& target) :
    _source(source),
    _target(target) {
}

inline Vector2 Segment::normal() const {
    Vector2 step = direction();
    double norm = step.norm();

    if (norm <= 0.0) {
        return Vector2::Zero();
    }

    return Vector2(step.y(), -step.x()) / norm;
}

inline Moments Segment::moments(int max_degree) const {
    Moments result = parametric_moments(max_degree);
    double scale = length();

    for (const MultiIndex& index : MultiIndex::all_up_to(max_degree)) {
        result.set(index, scale * result.at(index));
    }

    return result;
}

// Expanding (source + t * direction) binomially in each coordinate turns
// every monomial into a sum of powers of t, which integrate to 1 / (i + j + 1)
// over the unit parameter interval.
inline Moments Segment::parametric_moments(int max_degree) const {
    const std::vector<std::vector<double>>& binomial = globe::math::binomials(max_degree);
    Vector2 step = direction();
    Moments result(max_degree);

    std::vector<double> source_x = powers(_source.x(), max_degree);
    std::vector<double> source_y = powers(_source.y(), max_degree);
    std::vector<double> step_x = powers(step.x(), max_degree);
    std::vector<double> step_y = powers(step.y(), max_degree);

    for (const MultiIndex& index : MultiIndex::all_up_to(max_degree)) {
        if (index.z != 0) {
            continue;
        }

        double total = 0.0;

        for (int i = 0; i <= index.x; ++i) {
            for (int j = 0; j <= index.y; ++j) {
                total +=
                    binomial[index.x][i] * source_x[index.x - i] * step_x[i] *
                    binomial[index.y][j] * source_y[index.y - j] * step_y[j] /
                    static_cast<double>(i + j + 1);
            }
        }

        result.set(index, total);
    }

    return result;
}

} // namespace globe::geometry::planar

#endif //GLOBEART_SRC_GLOBE_GEOMETRY_PLANAR_SEGMENT_HPP_
