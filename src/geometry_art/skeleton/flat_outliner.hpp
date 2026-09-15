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
#include <variant>

namespace geometry_art::skeleton {

using geometry::planar::Embedding;
using geometry::planar::Segment;
using voronoi::flat::CellClipper;
using voronoi::flat::Diagram;
using voronoi::flat::Wall;

// Outlines the cells of a flat tessellation on the surface its domain
// embeds in. A cell's chart is the plane it is already computed in. Walls
// are pushed out by the half width so the bar along the domain's boundary
// is as wide as any other and centred on it; every other cut is inset.
class FlatOutliner {
 public:
    explicit FlatOutliner(const Diagram& diagram);

    [[nodiscard]] size_t size() const { return _diagram.size(); }
    [[nodiscard]] Outline outline(size_t index, double half_width, double step) const;

 private:
    const Diagram& _diagram;
    Embedding _embedding;

    [[nodiscard]] Loop sampled(const CellClipper& region, double step) const;
    [[nodiscard]] Sample sample(const Vector2& point, bool wall) const;
};

inline FlatOutliner::FlatOutliner(const Diagram& diagram) :
    _diagram(diagram),
    _embedding(Embedding::of(diagram.domain())) {
}

inline Outline FlatOutliner::outline(size_t index, double half_width, double step) const {
    return Outline{
        sampled(_diagram.offset_cell(index, 0.0, -half_width), step),
        sampled(_diagram.offset_cell(index, half_width, half_width), step)
    };
}

inline Loop FlatOutliner::sampled(const CellClipper& region, double step) const {
    Loop loop;

    for (size_t index = 0; index < region.size(); ++index) {
        const CellClipper::Edge& edge = region.edges()[index];
        Segment segment(edge.source, region.target(index));
        bool wall = std::holds_alternative<Wall>(edge.cut);
        size_t count = _embedding.curved() ? segments(segment.length(), step) : 1;

        for (size_t piece = 0; piece < count; ++piece) {
            loop.push_back(sample(segment.interpolate(static_cast<double>(piece) / count), wall));
        }
    }

    return loop;
}

inline Sample FlatOutliner::sample(const Vector2& point, bool wall) const {
    return Sample{point, _embedding.position(point), _embedding.normal(point), wall};
}

} // namespace geometry_art::skeleton

#endif //GEOMETRY_ART_SKELETON_FLAT_OUTLINER_HPP_
