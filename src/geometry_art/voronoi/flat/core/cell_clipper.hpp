#ifndef GEOMETRY_ART_VORONOI_FLAT_CORE_CELL_CLIPPER_HPP_
#define GEOMETRY_ART_VORONOI_FLAT_CORE_CELL_CLIPPER_HPP_

#include "cut.hpp"
#include "../../../types.hpp"
#include <cstddef>
#include <utility>
#include <vector>

namespace geometry_art::voronoi::flat {

// A convex cell under construction as an intersection of half planes, kept
// as a counter-clockwise loop of edges that each remember the cut they lie
// on. Clipping by one more half plane keeps every edge's cut and gives the
// closing edge the new one, so when the loop is finished each vertex knows
// the two cuts that meet there.
class CellClipper {
 public:
    struct Edge {
        Vector2 source;
        Cut cut;
    };

    explicit CellClipper(std::vector<Edge> edges);

    [[nodiscard]] const std::vector<Edge>& edges() const { return _edges; }
    [[nodiscard]] size_t size() const { return _edges.size(); }
    [[nodiscard]] const Vector2& target(size_t index) const { return _edges[(index + 1) % _edges.size()].source; }

    // Half planes are given as an inward normal and a point on the boundary,
    // matching how a Voronoi bisector is described.
    void clip(const Vector2& inward_normal, const Vector2& boundary_point, Cut cut);

 private:
    std::vector<Edge> _edges;

    [[nodiscard]] static Vector2 crossing(
        const Vector2& from,
        const Vector2& to,
        double from_offset,
        double to_offset
    );
};

inline CellClipper::CellClipper(std::vector<Edge> edges) :
    _edges(std::move(edges)) {
}

inline void CellClipper::clip(const Vector2& inward_normal, const Vector2& boundary_point, Cut cut) {
    std::vector<Edge> kept;
    kept.reserve(_edges.size() + 1);

    for (size_t index = 0; index < _edges.size(); ++index) {
        const Edge& edge = _edges[index];
        const Vector2& to = target(index);
        double from_offset = inward_normal.dot(edge.source - boundary_point);
        double to_offset = inward_normal.dot(to - boundary_point);
        bool from_inside = from_offset >= 0.0;
        bool to_inside = to_offset >= 0.0;

        if (from_inside) {
            kept.push_back(edge);
        }

        if (from_inside && !to_inside) {
            kept.push_back(Edge{crossing(edge.source, to, from_offset, to_offset), cut});
        }

        if (!from_inside && to_inside) {
            kept.push_back(Edge{crossing(edge.source, to, from_offset, to_offset), edge.cut});
        }
    }

    if (kept.size() < 3) {
        kept.clear();
    }

    _edges = std::move(kept);
}

inline Vector2 CellClipper::crossing(
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

} // namespace geometry_art::voronoi::flat

#endif //GEOMETRY_ART_VORONOI_FLAT_CORE_CELL_CLIPPER_HPP_
