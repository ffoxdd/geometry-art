#ifndef GEOMETRY_ART_GEOMETRY_SPHERICAL_POLYGON_HPP_
#define GEOMETRY_ART_GEOMETRY_SPHERICAL_POLYGON_HPP_

#include "../bounding_box.hpp"
#include "../arc.hpp"
#include "../helpers.hpp"
#include "bounding_box_calculator.hpp"
#include "../../../types.hpp"
#include "../../../math/polynomial/moments.hpp"
#include "../../../math/polynomial/multi_index.hpp"
#include "../../../std_ext/ranges.hpp"
#include <cstdint>
#include <Eigen/Core>
#include <utility>
#include <vector>
#include <cmath>
#include <ranges>
#include <algorithm>
#include <cassert>
#include <optional>

namespace geometry_art::geometry::spherical::polygon {

using geometry_art::VectorS2;
using geometry_art::std_ext::all_circular_adjacent_pairs;
using geometry_art::std_ext::circular_adjacent_pairs;
using geometry_art::geometry::spherical::Arc;
using geometry_art::geometry::spherical::BoundingBox;
using geometry_art::math::polynomial::Moments;
using geometry_art::math::polynomial::MultiIndex;

class Polygon {
 public:
    explicit Polygon(std::vector<Arc> arcs);

    [[nodiscard]] const std::vector<Arc>& arcs() const { return _arcs; }
    [[nodiscard]] auto points() const;
    [[nodiscard]] BoundingBox bounding_box() const;
    [[nodiscard]] VectorS2 centroid() const;
    [[nodiscard]] double bounding_sphere_radius() const;

    [[nodiscard]] double area() const;
    [[nodiscard]] Moments moments(int max_degree) const;
    [[nodiscard]] Moments moments(int max_degree, const std::vector<Moments>& arc_moments) const;
    [[nodiscard]] VectorS2 first_moment() const;
    [[nodiscard]] Eigen::Matrix3d second_moment() const;

    [[nodiscard]] bool contains(const VectorS2& point) const;
    [[nodiscard]] std::optional<Polygon> clipped_by(const VectorS2& half_space_normal) const;

    [[nodiscard]] static std::optional<Polygon> from_points(const std::vector<VectorS2>& points);

 private:
    std::vector<Arc> _arcs;

    [[nodiscard]] bool arcs_form_closed_loop() const;
    [[nodiscard]] std::vector<double> turning_angles() const;
    [[nodiscard]] VectorS2 vertex_average() const;

    [[nodiscard]] static double turning_angle(const Arc& incoming, const Arc& outgoing);
    [[nodiscard]] static int leading_axis(const MultiIndex& index);
};

inline Polygon::Polygon(std::vector<Arc> arcs) :
    _arcs(std::move(arcs)) {

    assert(!_arcs.empty());
    assert(arcs_form_closed_loop());
}

inline auto Polygon::points() const {
    return _arcs | std::views::transform(
        [](const Arc& arc) { return arc.source(); }
    );
}

inline bool Polygon::arcs_form_closed_loop() const {
    return all_circular_adjacent_pairs(_arcs, [](const auto& pair) {
        const auto& [previous, next] = pair;
        return (previous.target() - next.source()).squaredNorm() < GEOMETRIC_EPSILON;
    });
}

inline BoundingBox Polygon::bounding_box() const {
    return BoundingBoxCalculator(_arcs).calculate();
}

inline VectorS2 Polygon::centroid() const {
    VectorS2 moment = first_moment();
    double norm = moment.norm();

    if (norm < GEOMETRIC_EPSILON) {
        return vertex_average();
    }

    return moment / norm;
}

inline VectorS2 Polygon::vertex_average() const {
    VectorS2 sum = VectorS2::Zero();

    for (const VectorS2& point : points()) {
        sum += point;
    }

    double norm = sum.norm();
    return norm < GEOMETRIC_EPSILON ? bounding_box().center() : VectorS2(sum / norm);
}

inline double Polygon::bounding_sphere_radius() const {
    VectorS2 center = centroid();
    double max_squared_distance = 0.0;

    for (const VectorS2& point : points()) {
        max_squared_distance = std::max(max_squared_distance, (point - center).squaredNorm());
    }

    return std::sqrt(max_squared_distance);
}

inline double Polygon::area() const {
    std::vector<double> angles = turning_angles();
    std::sort(angles.begin(), angles.end());

    double turning_sum = 0.0;
    for (double angle : angles) {
        turning_sum += angle;
    }

    return TWO_PI - turning_sum;
}

inline std::vector<double> Polygon::turning_angles() const {
    std::vector<double> angles;
    angles.reserve(_arcs.size());

    for (const auto& pair : circular_adjacent_pairs(_arcs)) {
        const auto& [incoming, outgoing] = pair;
        angles.push_back(turning_angle(incoming, outgoing));
    }

    return angles;
}

inline double Polygon::turning_angle(const Arc& incoming, const Arc& outgoing) {
    const VectorS2& vertex = outgoing.source();
    VectorS2 direction_in = incoming.normal().cross(vertex);
    VectorS2 direction_out = outgoing.tangent_at_source();

    return std::atan2(direction_in.cross(direction_out).dot(vertex), direction_in.dot(direction_out));
}

inline Moments Polygon::moments(int max_degree) const {
    std::vector<Moments> arc_moments;
    arc_moments.reserve(_arcs.size());

    for (const Arc& arc : _arcs) {
        arc_moments.push_back(arc.moments(std::max(max_degree - 1, 0)));
    }

    return moments(max_degree, arc_moments);
}

inline Moments Polygon::moments(int max_degree, const std::vector<Moments>& arc_moments) const {
    assert(arc_moments.size() == _arcs.size());

    Moments result(max_degree);
    result.set(MultiIndex{0, 0, 0}, area());

    for (int degree = 1; degree <= max_degree; ++degree) {
        for (const MultiIndex& beta : MultiIndex::all_of_degree(degree)) {
            int axis = leading_axis(beta);
            MultiIndex alpha = beta.lowered(axis);

            double boundary_term = 0.0;
            for (size_t i = 0; i < _arcs.size(); ++i) {
                boundary_term += _arcs[i].normal()[axis] * arc_moments[i].at(alpha);
            }

            double interior_term = alpha[axis] > 0 ? alpha[axis] * result.at(alpha.lowered(axis)) : 0.0;

            result.set(beta, (interior_term + boundary_term) / (alpha.degree() + 2));
        }
    }

    return result;
}

inline int Polygon::leading_axis(const MultiIndex& index) {
    return index.x > 0 ? 0 : (index.y > 0 ? 1 : 2);
}

inline VectorS2 Polygon::first_moment() const {
    Moments moments = this->moments(1);
    return VectorS2(moments.at(1, 0, 0), moments.at(0, 1, 0), moments.at(0, 0, 1));
}

inline Eigen::Matrix3d Polygon::second_moment() const {
    Moments moments = this->moments(2);
    Eigen::Matrix3d result;

    for (int row = 0; row < 3; ++row) {
        for (int column = 0; column < 3; ++column) {
            result(row, column) = moments.at(MultiIndex::unit(row).raised(column));
        }
    }

    return result;
}

inline bool Polygon::contains(const VectorS2& point) const {
    for (const auto& arc : _arcs) {
        if (arc.normal().dot(point) < -GEOMETRIC_EPSILON) {
            return false;
        }
    }

    return true;
}

// Sutherland-Hodgman against one hemisphere. Every surviving arc keeps the
// great circle it already lies on, and the arcs closing the polygon along
// the clip circle take that circle's normal, so no normal is ever rederived
// from two nearby points.
inline std::optional<Polygon> Polygon::clipped_by(const VectorS2& half_space_normal) const {
    std::vector<Arc> kept;
    std::optional<VectorS2> last_exit;
    std::optional<VectorS2> first_entry;
    bool all_inside = true;

    for (const Arc& arc : _arcs) {
        bool source_inside = half_space_normal.dot(arc.source()) >= -GEOMETRIC_EPSILON;
        bool target_inside = half_space_normal.dot(arc.target()) >= -GEOMETRIC_EPSILON;
        all_inside = all_inside && source_inside && target_inside;

        if (source_inside && target_inside) {
            kept.push_back(arc);
            continue;
        }

        if (!source_inside && !target_inside) {
            continue;
        }

        VectorS2 crossing = arc.crossing_with(half_space_normal);

        if (source_inside) {
            kept.emplace_back(arc.source(), crossing, arc.normal());
            last_exit = crossing;
            continue;
        }

        if (last_exit) {
            kept.emplace_back(*last_exit, crossing, half_space_normal);
        } else {
            first_entry = crossing;
        }

        kept.emplace_back(crossing, arc.target(), arc.normal());
    }

    if (all_inside) {
        return *this;
    }

    if (last_exit && first_entry) {
        kept.emplace_back(*last_exit, *first_entry, half_space_normal);
    }

    if (kept.size() < 2) {
        return std::nullopt;
    }

    return Polygon(std::move(kept));
}

inline std::optional<Polygon> Polygon::from_points(const std::vector<VectorS2>& points) {
    std::vector<VectorS2> distinct;

    for (const VectorS2& point : points) {
        if (distinct.empty() || (point - distinct.back()).squaredNorm() > GEOMETRIC_EPSILON * GEOMETRIC_EPSILON) {
            distinct.push_back(point);
        }
    }

    while (distinct.size() > 1 && (distinct.front() - distinct.back()).squaredNorm() <= GEOMETRIC_EPSILON * GEOMETRIC_EPSILON) {
        distinct.pop_back();
    }

    if (distinct.size() < 3) {
        return std::nullopt;
    }

    std::vector<Arc> arcs;
    arcs.reserve(distinct.size());
    for (size_t i = 0; i < distinct.size(); ++i) {
        arcs.emplace_back(distinct[i], distinct[(i + 1) % distinct.size()]);
    }

    return Polygon(std::move(arcs));
}

} // namespace geometry_art::geometry::spherical::polygon

namespace geometry_art::geometry::spherical {
using Polygon = polygon::Polygon;
}

namespace geometry_art {
using Polygon = geometry::spherical::polygon::Polygon;
}

#endif //GEOMETRY_ART_GEOMETRY_SPHERICAL_POLYGON_HPP_
