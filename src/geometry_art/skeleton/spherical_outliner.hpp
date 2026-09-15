#ifndef GEOMETRY_ART_SKELETON_SPHERICAL_OUTLINER_HPP_
#define GEOMETRY_ART_SKELETON_SPHERICAL_OUTLINER_HPP_

#include "outline.hpp"
#include "../geometry/spherical/arc.hpp"
#include "../geometry/spherical/cap_polygon.hpp"
#include "../geometry/spherical/polygon/polygon.hpp"
#include "../types.hpp"
#include "../voronoi/spherical/core/sphere.hpp"
#include <cmath>
#include <cstddef>

namespace geometry_art::skeleton {

using geometry::spherical::Arc;
using geometry::spherical::CapPolygon;
using geometry::spherical::Polygon;
using voronoi::spherical::Sphere;

// Outlines the cells of a spherical tessellation. The sphere is its own
// embedding, with the normal the position itself, and each cell is charted
// gnomonically about its centroid, which keeps its geodesic edges straight
// and its inset convex.
class SphericalOutliner {
 public:
    explicit SphericalOutliner(const Sphere& sphere);

    [[nodiscard]] size_t size() const { return _sphere.cell_count(); }
    [[nodiscard]] Outline outline(size_t index, double half_width, double step) const;

 private:
    struct Chart {
        VectorS2 center;
        VectorS2 first;
        VectorS2 second;

        [[nodiscard]] Sample sample(const VectorS2& point) const;
    };

    const Sphere& _sphere;

    [[nodiscard]] static Chart chart_about(const VectorS2& center);
    [[nodiscard]] static VectorS2 perpendicular_to(const VectorS2& axis);
};

inline SphericalOutliner::SphericalOutliner(const Sphere& sphere) :
    _sphere(sphere) {
}

inline Outline SphericalOutliner::outline(size_t index, double half_width, double step) const {
    Polygon cell = _sphere.cell(index);
    CapPolygon inset = _sphere.offset_cell(index, half_width);
    Chart chart = chart_about(cell.centroid());
    Outline outline;

    for (const Arc& arc : cell.arcs()) {
        size_t count = segments(arc.length(), step);

        for (size_t piece = 0; piece < count; ++piece) {
            outline.outer.push_back(chart.sample(arc.interpolate(static_cast<double>(piece) / count)));
        }
    }

    for (size_t edge = 0; edge < inset.size(); ++edge) {
        size_t count = segments(inset.length(edge), step);

        for (size_t piece = 0; piece < count; ++piece) {
            outline.inner.push_back(chart.sample(inset.point(edge, static_cast<double>(piece) / count)));
        }
    }

    return outline;
}

inline Sample SphericalOutliner::Chart::sample(const VectorS2& point) const {
    double depth = point.dot(center);
    return Sample{Vector2(point.dot(first) / depth, point.dot(second) / depth), point, point, false};
}

inline SphericalOutliner::Chart SphericalOutliner::chart_about(const VectorS2& center) {
    VectorS2 first = perpendicular_to(center);
    return Chart{center, first, center.cross(first)};
}

inline VectorS2 SphericalOutliner::perpendicular_to(const VectorS2& axis) {
    VectorS2 candidate = (std::abs(axis.z()) < 0.9) ? VectorS2(0.0, 0.0, 1.0) : VectorS2(1.0, 0.0, 0.0);
    return (candidate - axis.dot(candidate) * axis).normalized();
}

} // namespace geometry_art::skeleton

#endif //GEOMETRY_ART_SKELETON_SPHERICAL_OUTLINER_HPP_
