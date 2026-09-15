#ifndef GEOMETRY_ART_GEOMETRY_SPHERICAL_CAP_POLYGON_HPP_
#define GEOMETRY_ART_GEOMETRY_SPHERICAL_CAP_POLYGON_HPP_

#include "cap.hpp"
#include "polygon/polygon.hpp"
#include "../../types.hpp"
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

namespace geometry_art::geometry::spherical {

// An intersection of caps, kept as a loop of edges that each run along the
// rim of the cap they came from, counter-clockwise about its axis with the
// region on the left. A great-circle polygon is the case where every cap is
// a hemisphere; shrinking one by a distance replaces its great circles with
// small ones, which is why the edges are rim arcs rather than geodesics.
//
// Clipping by one more cap keeps the parts of every edge inside it. A rim
// arc is not a geodesic, so the parts inside are found on the rim itself:
// along a rim the new cap covers one angular interval, and an edge keeps
// whatever of its own interval falls in it, which may be two pieces.
class CapPolygon {
 public:
    struct Edge {
        VectorS2 source;
        Cap cap;
    };

    explicit CapPolygon(std::vector<Edge> edges);

    [[nodiscard]] static CapPolygon of(const polygon::Polygon& polygon);

    [[nodiscard]] const std::vector<Edge>& edges() const { return _edges; }
    [[nodiscard]] size_t size() const { return _edges.size(); }
    [[nodiscard]] bool empty() const { return _edges.empty(); }
    [[nodiscard]] const VectorS2& target(size_t index) const { return _edges[(index + 1) % _edges.size()].source; }

    [[nodiscard]] double length(size_t index) const;
    [[nodiscard]] VectorS2 point(size_t index, double fraction) const;

    void clip(const Cap& cap);

 private:
    static constexpr double ANGLE_EPSILON = 1e-9;

    struct Rim {
        VectorS2 center;
        VectorS2 first;
        VectorS2 second;
        double radius;

        [[nodiscard]] VectorS2 point(double angle) const;
        [[nodiscard]] double angle_of(const VectorS2& point) const;
    };

    struct Interval {
        double start;
        double end;
    };

    std::vector<Edge> _edges;

    [[nodiscard]] double turn(size_t index) const;
    [[nodiscard]] std::vector<Interval> inside(size_t index, const Cap& cap) const;

    [[nodiscard]] static Rim rim_of(const Cap& cap);
    [[nodiscard]] static VectorS2 perpendicular_to(const VectorS2& axis);
    [[nodiscard]] static double wrapped(double angle);
};

inline CapPolygon::CapPolygon(std::vector<Edge> edges) :
    _edges(std::move(edges)) {
}

inline CapPolygon CapPolygon::of(const polygon::Polygon& polygon) {
    std::vector<Edge> edges;
    edges.reserve(polygon.arcs().size());

    for (const Arc& arc : polygon.arcs()) {
        edges.push_back(Edge{arc.source(), Cap::hemisphere(arc.normal())});
    }

    return CapPolygon(std::move(edges));
}

inline double CapPolygon::length(size_t index) const {
    return turn(index) * rim_of(_edges[index].cap).radius;
}

inline VectorS2 CapPolygon::point(size_t index, double fraction) const {
    Rim rim = rim_of(_edges[index].cap);
    return rim.point(rim.angle_of(_edges[index].source) + fraction * turn(index));
}

inline void CapPolygon::clip(const Cap& cap) {
    std::vector<Edge> kept;
    kept.reserve(_edges.size() + 2);

    for (size_t index = 0; index < _edges.size(); ++index) {
        const Edge& edge = _edges[index];
        Rim rim = rim_of(edge.cap);
        double start = rim.angle_of(edge.source);
        double turn = this->turn(index);

        for (const Interval& piece : inside(index, cap)) {
            if (piece.start <= ANGLE_EPSILON) {
                kept.push_back(edge);
            } else {
                kept.push_back(Edge{rim.point(start + piece.start), edge.cap});
            }

            if (piece.end < turn - ANGLE_EPSILON) {
                kept.push_back(Edge{rim.point(start + piece.end), cap});
            }
        }
    }

    if (kept.size() < 2) {
        kept.clear();
    }

    _edges = std::move(kept);
}

// The angle the edge sweeps about its axis. A lone edge is a whole rim.
inline double CapPolygon::turn(size_t index) const {
    Rim rim = rim_of(_edges[index].cap);
    double sweep = wrapped(rim.angle_of(target(index)) - rim.angle_of(_edges[index].source));

    if (sweep < ANGLE_EPSILON && _edges.size() == 1) {
        return TWO_PI;
    }

    return sweep;
}

// Along the edge's rim the cap's height is a sinusoid in the rim angle, so
// the rim is inside the cap on one angular interval, centred where the
// sinusoid peaks. The pieces of the edge in that interval are reported as
// angles from the edge's source, in travel order.
inline std::vector<CapPolygon::Interval> CapPolygon::inside(size_t index, const Cap& cap) const {
    const Edge& edge = _edges[index];
    Rim rim = rim_of(edge.cap);
    double sweep = turn(index);
    double base = cap.axis.dot(rim.center) - cap.offset;
    double along_first = rim.radius * cap.axis.dot(rim.first);
    double along_second = rim.radius * cap.axis.dot(rim.second);
    double amplitude = std::hypot(along_first, along_second);

    if (amplitude < GEOMETRIC_EPSILON) {
        return base >= 0.0 ? std::vector<Interval>{Interval{0.0, sweep}} : std::vector<Interval>{};
    }

    double threshold = -base / amplitude;

    if (threshold <= -1.0) {
        return {Interval{0.0, sweep}};
    }

    if (threshold >= 1.0) {
        return {};
    }

    double half_width = std::acos(threshold);
    double peak = std::atan2(along_second, along_first);
    double entry = wrapped(peak - half_width - rim.angle_of(edge.source));
    std::vector<Interval> pieces;

    for (double start : {entry - TWO_PI, entry}) {
        Interval piece{std::max(start, 0.0), std::min(start + 2.0 * half_width, sweep)};

        if (piece.end - piece.start > ANGLE_EPSILON) {
            pieces.push_back(piece);
        }
    }

    return pieces;
}

inline VectorS2 CapPolygon::Rim::point(double angle) const {
    return center + radius * (std::cos(angle) * first + std::sin(angle) * second);
}

inline double CapPolygon::Rim::angle_of(const VectorS2& point) const {
    return std::atan2(point.dot(second), point.dot(first));
}

inline CapPolygon::Rim CapPolygon::rim_of(const Cap& cap) {
    VectorS2 first = perpendicular_to(cap.axis);
    return Rim{cap.rim_center(), first, cap.axis.cross(first), cap.rim_radius()};
}

inline VectorS2 CapPolygon::perpendicular_to(const VectorS2& axis) {
    VectorS2 candidate = (std::abs(axis.z()) < 0.9) ? VectorS2(0.0, 0.0, 1.0) : VectorS2(1.0, 0.0, 0.0);
    return (candidate - axis.dot(candidate) * axis).normalized();
}

inline double CapPolygon::wrapped(double angle) {
    double reduced = std::fmod(angle, TWO_PI);
    return reduced < 0.0 ? reduced + TWO_PI : reduced;
}

} // namespace geometry_art::geometry::spherical

#endif //GEOMETRY_ART_GEOMETRY_SPHERICAL_CAP_POLYGON_HPP_
