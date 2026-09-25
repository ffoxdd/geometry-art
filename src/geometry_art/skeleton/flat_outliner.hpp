#ifndef GEOMETRY_ART_SKELETON_FLAT_OUTLINER_HPP_
#define GEOMETRY_ART_SKELETON_FLAT_OUTLINER_HPP_

#include "outline.hpp"
#include "../geometry/planar/embedding.hpp"
#include "../geometry/planar/segment.hpp"
#include "../types.hpp"
#include "../voronoi/flat/core/cell_clipper.hpp"
#include "../voronoi/flat/core/cut.hpp"
#include "../voronoi/flat/core/diagram.hpp"
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <utility>
#include <variant>
#include <vector>

namespace geometry_art::skeleton {

using geometry::planar::Domain;
using geometry::planar::Embedding;
using geometry::planar::Segment;
using voronoi::flat::CellClipper;
using voronoi::flat::Diagram;
using voronoi::flat::Wall;

// A rectangle centred on the domain that the skeleton is cut down to, laid
// flat, with a frame of the given width running around it outside. What
// shows through the window is the tessellation away from the domain's own
// boundary.
struct Window {
    Vector2 size;
    double frame_width;
};

// Outlines the cells of a flat tessellation on the surface its domain
// embeds in. A cell's chart is the plane it is already computed in. Walls
// are pushed out by the half width so the bar along the domain's boundary
// is as wide as any other and centred on it; every other cut is inset.
// Through a window, each cell is clipped to it with its bars' openings
// ending at the window's edge and its ring running on to the frame's
// outer edge, which is walled. Along a wrapped axis the window recurs every
// period and a cell's chart can reach the copies on either side of it, so
// a cell is outlined once per copy, and a copy it misses outlines nothing.
class FlatOutliner {
 public:
    explicit FlatOutliner(const Diagram& diagram, std::optional<Window> window = std::nullopt);

    [[nodiscard]] size_t size() const { return _diagram.size() * _copies.size(); }
    [[nodiscard]] Outline outline(size_t index, double half_width, double step) const;

 private:
    const Diagram& _diagram;
    std::optional<Window> _window;
    std::vector<Vector2> _copies;
    Embedding _embedding;

    [[nodiscard]] Loop sampled(const CellClipper& region, double step, const Vector2& shift) const;
    [[nodiscard]] Sample sample(const Vector2& point, bool wall) const;

    static void clip_to_rectangle(CellClipper& region, const Vector2& center, const Vector2& half_size);
    [[nodiscard]] static std::optional<Window> fitted(const std::optional<Window>& window, const Domain& domain);
    [[nodiscard]] static std::vector<Vector2> copies_of(const Domain& domain, const std::optional<Window>& window);
    [[nodiscard]] static Embedding embedding_of(const Domain& domain, const std::optional<Window>& window);
};

inline FlatOutliner::FlatOutliner(const Diagram& diagram, std::optional<Window> window) :
    _diagram(diagram),
    _window(fitted(window, diagram.domain())),
    _copies(copies_of(diagram.domain(), window)),
    _embedding(embedding_of(diagram.domain(), window)) {
}

inline Outline FlatOutliner::outline(size_t index, double half_width, double step) const {
    size_t cell = index / _copies.size();
    const Vector2& shift = _copies[index % _copies.size()];
    CellClipper outer = _diagram.offset_cell(cell, 0.0, -half_width);
    CellClipper inner = _diagram.offset_cell(cell, half_width, half_width);

    if (!_window) {
        return Outline{sampled(outer, step, shift), sampled(inner, step, shift)};
    }

    Vector2 center = 0.5 * Vector2(_diagram.domain().width, _diagram.domain().height) + shift;
    Vector2 half_size = 0.5 * _window->size;

    clip_to_rectangle(outer, center, half_size + Vector2::Constant(_window->frame_width));
    clip_to_rectangle(inner, center, half_size);

    return Outline{sampled(outer, step, shift), sampled(inner, step, shift)};
}

inline void FlatOutliner::clip_to_rectangle(CellClipper& region, const Vector2& center, const Vector2& half_size) {
    for (int axis = 0; axis < 2; ++axis) {
        for (double side : {-1.0, 1.0}) {
            Vector2 inward = Vector2::Zero();
            inward[axis] = -side;
            Vector2 boundary = center;
            boundary[axis] += side * half_size[axis];
            region.clip(inward, boundary, Wall{inward});
        }
    }
}

inline std::optional<Window> FlatOutliner::fitted(const std::optional<Window>& window, const Domain& domain) {
    if (!window) {
        return window;
    }

    Vector2 framed = window->size + Vector2::Constant(2.0 * window->frame_width);

    if (framed.x() > domain.width || framed.y() > domain.height) {
        throw std::invalid_argument("the window and its frame must fit within the domain");
    }

    return window;
}

// A cell's chart stays within half a period of its site, so along a wrapped
// axis it can reach the window's copies one period either side.
inline std::vector<Vector2> FlatOutliner::copies_of(const Domain& domain, const std::optional<Window>& window) {
    std::vector<Vector2> copies{Vector2::Zero()};

    if (!window) {
        return copies;
    }

    for (int axis = 0; axis < 2; ++axis) {
        if (!domain.wrapped(axis)) {
            continue;
        }

        std::vector<Vector2> spread;

        for (const Vector2& copy : copies) {
            for (double period : {-1.0, 0.0, 1.0}) {
                Vector2 shifted = copy;
                shifted[axis] += period * domain.extent(axis);
                spread.push_back(shifted);
            }
        }

        copies = std::move(spread);
    }

    return copies;
}

inline Embedding FlatOutliner::embedding_of(const Domain& domain, const std::optional<Window>& window) {
    if (!window) {
        return Embedding::of(domain);
    }

    return Embedding::of(Domain::plane(domain.width, domain.height));
}

inline Loop FlatOutliner::sampled(const CellClipper& region, double step, const Vector2& shift) const {
    Loop loop;

    for (size_t index = 0; index < region.size(); ++index) {
        const CellClipper::Edge& edge = region.edges()[index];
        Segment segment(edge.source, region.target(index));
        bool wall = std::holds_alternative<Wall>(edge.cut);
        size_t count = _embedding.curved() ? segments(segment.length(), step) : 1;

        for (size_t piece = 0; piece < count; ++piece) {
            loop.push_back(sample(segment.interpolate(static_cast<double>(piece) / count) - shift, wall));
        }
    }

    return loop;
}

inline Sample FlatOutliner::sample(const Vector2& point, bool wall) const {
    return Sample{point, _embedding.position(point), _embedding.normal(point), wall};
}

} // namespace geometry_art::skeleton

#endif //GEOMETRY_ART_SKELETON_FLAT_OUTLINER_HPP_
