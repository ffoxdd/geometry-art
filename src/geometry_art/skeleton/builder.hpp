#ifndef GEOMETRY_ART_SKELETON_BUILDER_HPP_
#define GEOMETRY_ART_SKELETON_BUILDER_HPP_

#include "annulus.hpp"
#include "outline.hpp"
#include "welder.hpp"
#include "../io/mesh/types.hpp"
#include "../types.hpp"
#include <cstddef>
#include <vector>

namespace geometry_art::skeleton {

// The bars in the model's own units: how wide they run along the surface,
// how thick they stand off it, and how finely curved edges are sampled.
// The scale multiplies every position on the way out.
struct Parameters {
    double bar_width;
    double bar_thickness;
    double resolution;
    double scale = 1.0;
};

// Thickens the edge graph of a tessellation into a solid. The region within
// half a bar of the graph is the domain less every cell's inset, so each
// cell contributes the ring between its boundary and its inset, and the
// rings tile the surface with mitred joints wherever cells meet. The ring
// is laid at both offsets along the normal, its inset walled around, and
// its boundary walled wherever it is the domain's own.
class Builder {
 public:
    explicit Builder(Parameters parameters);

    template<Outliner OutlinerType>
    [[nodiscard]] SurfaceMesh build(const OutlinerType& outliner) const;

 private:
    Parameters _parameters;

    void add_cell(const Outline& outline, Welder& welder) const;
    void add_wall(const Sample& from, const Sample& to, Welder& welder) const;

    [[nodiscard]] Vector3 top(const Sample& sample) const;
    [[nodiscard]] Vector3 bottom(const Sample& sample) const;
};

inline Builder::Builder(Parameters parameters) :
    _parameters(parameters) {
}

template<Outliner OutlinerType>
SurfaceMesh Builder::build(const OutlinerType& outliner) const {
    Welder welder(_parameters.scale);

    for (size_t index = 0; index < outliner.size(); ++index) {
        add_cell(outliner.outline(index, 0.5 * _parameters.bar_width, _parameters.resolution), welder);
    }

    return welder.take();
}

inline void Builder::add_cell(const Outline& outline, Welder& welder) const {
    if (outline.outer.size() < 3) {
        return;
    }

    std::vector<Sample> samples = outline.outer;
    samples.insert(samples.end(), outline.inner.begin(), outline.inner.end());

    std::vector<Vector2> outer_chart;
    std::vector<Vector2> inner_chart;

    for (const Sample& sample : outline.outer) {
        outer_chart.push_back(sample.chart);
    }

    for (const Sample& sample : outline.inner) {
        inner_chart.push_back(sample.chart);
    }

    std::vector<VertexIndex> tops;
    std::vector<VertexIndex> bottoms;

    for (const Sample& sample : samples) {
        tops.push_back(welder.vertex(top(sample)));
        bottoms.push_back(welder.vertex(bottom(sample)));
    }

    for (const Triangle& triangle : triangulate_annulus(outer_chart, inner_chart)) {
        welder.triangle(tops[triangle[0]], tops[triangle[1]], tops[triangle[2]]);
        welder.triangle(bottoms[triangle[0]], bottoms[triangle[2]], bottoms[triangle[1]]);
    }

    for (size_t index = 0; index < outline.inner.size(); ++index) {
        add_wall(outline.inner[index], outline.inner[(index + 1) % outline.inner.size()], welder);
    }

    for (size_t index = 0; index < outline.outer.size(); ++index) {
        if (outline.outer[index].wall) {
            add_wall(outline.outer[(index + 1) % outline.outer.size()], outline.outer[index], welder);
        }
    }
}

// A wall stands along an edge with the void on its left, facing it.
inline void Builder::add_wall(const Sample& from, const Sample& to, Welder& welder) const {
    VertexIndex top_from = welder.vertex(top(from));
    VertexIndex top_to = welder.vertex(top(to));
    VertexIndex bottom_from = welder.vertex(bottom(from));
    VertexIndex bottom_to = welder.vertex(bottom(to));

    welder.triangle(top_from, top_to, bottom_to);
    welder.triangle(top_from, bottom_to, bottom_from);
}

inline Vector3 Builder::top(const Sample& sample) const {
    return sample.position + 0.5 * _parameters.bar_thickness * sample.normal;
}

inline Vector3 Builder::bottom(const Sample& sample) const {
    return sample.position - 0.5 * _parameters.bar_thickness * sample.normal;
}

} // namespace geometry_art::skeleton

#endif //GEOMETRY_ART_SKELETON_BUILDER_HPP_
