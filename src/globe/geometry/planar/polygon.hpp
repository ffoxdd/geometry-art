#ifndef GLOBEART_SRC_GLOBE_GEOMETRY_PLANAR_POLYGON_HPP_
#define GLOBEART_SRC_GLOBE_GEOMETRY_PLANAR_POLYGON_HPP_

#include "segment.hpp"
#include "../../types.hpp"
#include "../../math/polynomial/moments.hpp"
#include "../../math/polynomial/multi_index.hpp"
#include <cstddef>
#include <optional>
#include <utility>
#include <vector>

namespace globe::geometry::planar {

// A closed loop of straight edges, oriented counter-clockwise so that the
// interior lies to the left. Moments follow the divergence theorem: an area
// integral of a monomial becomes a boundary integral of the monomial raised
// in x, which each edge integrates in closed form.
class Polygon {
 public:
    explicit Polygon(std::vector<Vector2> vertices);

    [[nodiscard]] const std::vector<Vector2>& vertices() const { return _vertices; }
    [[nodiscard]] size_t size() const { return _vertices.size(); }
    [[nodiscard]] std::vector<Segment> edges() const;

    [[nodiscard]] double area() const;
    [[nodiscard]] Vector2 centroid() const;
    [[nodiscard]] Moments moments(int max_degree) const;
    [[nodiscard]] bool contains(const Vector2& point) const;

    // Half planes are given as an inward normal and a point on the boundary,
    // matching how a Voronoi bisector is described.
    [[nodiscard]] std::optional<Polygon> clipped_by(const Vector2& inward_normal, const Vector2& boundary_point) const;

    [[nodiscard]] static Polygon rectangle(const Vector2& low, const Vector2& high);

 private:
    std::vector<Vector2> _vertices;

    [[nodiscard]] static Vector2 crossing(
        const Vector2& from,
        const Vector2& to,
        double from_offset,
        double to_offset
    );
};

inline Polygon::Polygon(std::vector<Vector2> vertices) :
    _vertices(std::move(vertices)) {
}

inline std::vector<Segment> Polygon::edges() const {
    std::vector<Segment> result;
    result.reserve(_vertices.size());

    for (size_t index = 0; index < _vertices.size(); ++index) {
        result.emplace_back(_vertices[index], _vertices[(index + 1) % _vertices.size()]);
    }

    return result;
}

inline double Polygon::area() const {
    double twice_area = 0.0;

    for (const Segment& edge : edges()) {
        twice_area += edge.source().x() * edge.target().y() - edge.target().x() * edge.source().y();
    }

    return 0.5 * twice_area;
}

inline Vector2 Polygon::centroid() const {
    Moments first = moments(1);
    double mass = first.at(0, 0, 0);

    if (mass == 0.0) {
        return Vector2::Zero();
    }

    return Vector2(first.at(1, 0, 0), first.at(0, 1, 0)) / mass;
}

// The divergence theorem with the field (x^(a+1) y^b / (a + 1), 0): the
// outward normal's x component times arc length is the edge direction's y
// component times the parameter step, so the length factors cancel.
inline Moments Polygon::moments(int max_degree) const {
    std::vector<Moments> edge_moments;

    for (const Segment& edge : edges()) {
        edge_moments.push_back(edge.parametric_moments(max_degree + 1));
    }

    Moments result(max_degree);

    for (const MultiIndex& index : MultiIndex::all_up_to(max_degree)) {
        if (index.z != 0) {
            continue;
        }

        MultiIndex raised = index.raised(0);
        double total = 0.0;
        size_t position = 0;

        for (const Segment& edge : edges()) {
            total += edge.direction().y() * edge_moments[position].at(raised);
            ++position;
        }

        result.set(index, total / static_cast<double>(index.x + 1));
    }

    return result;
}

inline bool Polygon::contains(const Vector2& point) const {
    for (const Segment& edge : edges()) {
        Vector2 step = edge.direction();
        Vector2 offset = point - edge.source();

        if (step.x() * offset.y() - step.y() * offset.x() < 0.0) {
            return false;
        }
    }

    return true;
}

inline std::optional<Polygon> Polygon::clipped_by(
    const Vector2& inward_normal,
    const Vector2& boundary_point
) const {
    std::vector<Vector2> kept;
    kept.reserve(_vertices.size() + 1);

    for (size_t index = 0; index < _vertices.size(); ++index) {
        const Vector2& from = _vertices[index];
        const Vector2& to = _vertices[(index + 1) % _vertices.size()];
        double from_offset = inward_normal.dot(from - boundary_point);
        double to_offset = inward_normal.dot(to - boundary_point);

        if (from_offset >= 0.0) {
            kept.push_back(from);
        }

        if ((from_offset >= 0.0) != (to_offset >= 0.0)) {
            kept.push_back(crossing(from, to, from_offset, to_offset));
        }
    }

    if (kept.size() < 3) {
        return std::nullopt;
    }

    return Polygon(std::move(kept));
}

inline Polygon Polygon::rectangle(const Vector2& low, const Vector2& high) {
    return Polygon(std::vector<Vector2>{
        Vector2(low.x(), low.y()),
        Vector2(high.x(), low.y()),
        Vector2(high.x(), high.y()),
        Vector2(low.x(), high.y())
    });
}

inline Vector2 Polygon::crossing(
    const Vector2& from,
    const Vector2& to,
    double from_offset,
    double to_offset
) {
    double span = from_offset - to_offset;

    if (span == 0.0) {
        return from;
    }

    return from + (from_offset / span) * (to - from);
}

} // namespace globe::geometry::planar

#endif //GLOBEART_SRC_GLOBE_GEOMETRY_PLANAR_POLYGON_HPP_
