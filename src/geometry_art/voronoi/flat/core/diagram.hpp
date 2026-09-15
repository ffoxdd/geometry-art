#ifndef GEOMETRY_ART_VORONOI_FLAT_CORE_DIAGRAM_HPP_
#define GEOMETRY_ART_VORONOI_FLAT_CORE_DIAGRAM_HPP_

#include "cell_clipper.hpp"
#include "cut.hpp"
#include "../../../types.hpp"
#include "../../../geometry/planar/domain.hpp"
#include "../../../geometry/planar/domain_delaunay.hpp"
#include "../../../geometry/planar/polygon.hpp"
#include "../../../geometry/planar/segment.hpp"
#include <cstddef>
#include <memory>
#include <utility>
#include <variant>
#include <vector>

namespace geometry_art::voronoi::flat {

using geometry::planar::Domain;
using geometry::planar::DomainDelaunay;
using geometry::planar::Neighbor;
using geometry::planar::Polygon;
using geometry::planar::Segment;

// One bisector edge of a cell. Everything geometric is expressed in the
// cell's own chart -- the plane around the cell's canonical site -- so the
// neighbour's position is carried here rather than read from its own
// canonical placement, which may sit a period away. The cuts at either end
// say what else pins that endpoint: a third site, or a wall.
struct CellEdgeInfo {
    size_t neighbor_index;
    Vector2 neighbor_position;
    Cut source_cut;
    Cut target_cut;
    Segment boundary;
};

// A Voronoi diagram on a flat domain, read off a Delaunay triangulation
// that is periodic along the wrapped axes. A cell is the site's fundamental
// region -- bounded by the bisectors with its own period images along
// wrapped axes and by the walls along walled ones -- clipped by the
// bisectors of its Delaunay neighbours. Bisector edges carry the sweep the
// optimizer differentiates; wall edges and self-image seams never move
// with the sites and only shape the cell.
class Diagram {
 public:
    Diagram(Domain domain, std::vector<Vector2> sites);

    Diagram(const Diagram&) = delete;
    Diagram& operator=(const Diagram&) = delete;
    Diagram(Diagram&&) = default;
    Diagram& operator=(Diagram&&) = default;

    [[nodiscard]] size_t size() const { return _triangulation.size(); }
    [[nodiscard]] const Domain& domain() const { return _triangulation.domain(); }
    [[nodiscard]] double width() const { return domain().width; }
    [[nodiscard]] double height() const { return domain().height; }
    [[nodiscard]] double area() const { return domain().area(); }

    [[nodiscard]] Vector2 site(size_t index) const { return _triangulation.site(index); }
    [[nodiscard]] Vector3 site_vector(size_t index) const;
    [[nodiscard]] std::vector<CellEdgeInfo> cell_edges(size_t index) const;
    [[nodiscard]] Polygon cell(size_t index) const;

    // The cell with each cut that bounds it moved inward by a distance
    // chosen per kind: bisectors and seams by one, walls by the other. A
    // negative distance moves the cut outward, and a cell its cuts consume
    // comes back empty.
    [[nodiscard]] CellClipper offset_cell(size_t index, double bisector_inset, double wall_inset) const;

    [[nodiscard]] Vector2 canonical(const Vector2& point) const { return domain().canonical(point); }

    [[nodiscard]] std::unique_ptr<Diagram> rebuilt(const std::vector<Vector3>& points) const;

 private:
    struct Inset {
        double bisector;
        double wall;
    };

    DomainDelaunay _triangulation;

    [[nodiscard]] CellClipper clipped_cell(size_t index, Inset inset) const;
    [[nodiscard]] CellClipper fundamental_region(size_t index, Inset inset) const;
};

inline Diagram::Diagram(Domain domain, std::vector<Vector2> sites) :
    _triangulation(domain, std::move(sites)) {
}

inline Vector3 Diagram::site_vector(size_t index) const {
    Vector2 point = site(index);
    return Vector3(point.x(), point.y(), 0.0);
}

inline std::vector<CellEdgeInfo> Diagram::cell_edges(size_t index) const {
    CellClipper cell = clipped_cell(index, Inset{0.0, 0.0});
    const std::vector<CellClipper::Edge>& edges = cell.edges();
    std::vector<CellEdgeInfo> result;
    result.reserve(edges.size());

    size_t count = edges.size();

    for (size_t position = 0; position < count; ++position) {
        const Bisector* bisector = std::get_if<Bisector>(&edges[position].cut);

        if (bisector == nullptr) {
            continue;
        }

        result.push_back(CellEdgeInfo{
            bisector->neighbor_index,
            bisector->neighbor_position,
            edges[(position + count - 1) % count].cut,
            edges[(position + 1) % count].cut,
            Segment(edges[position].source, cell.target(position))
        });
    }

    return result;
}

inline Polygon Diagram::cell(size_t index) const {
    CellClipper cell = clipped_cell(index, Inset{0.0, 0.0});
    std::vector<Vector2> vertices;
    vertices.reserve(cell.size());

    for (const CellClipper::Edge& edge : cell.edges()) {
        vertices.push_back(edge.source);
    }

    return Polygon(std::move(vertices));
}

inline CellClipper Diagram::offset_cell(size_t index, double bisector_inset, double wall_inset) const {
    return clipped_cell(index, Inset{bisector_inset, wall_inset});
}

inline std::unique_ptr<Diagram> Diagram::rebuilt(const std::vector<Vector3>& points) const {
    std::vector<Vector2> sites;
    sites.reserve(points.size());

    for (const Vector3& point : points) {
        sites.emplace_back(point.x(), point.y());
    }

    return std::make_unique<Diagram>(domain(), std::move(sites));
}

// Period images of the site itself are skipped: their bisectors are the
// fundamental region's own sides along wrapped axes, and the diagonal
// images' bisectors only touch its corners.
inline CellClipper Diagram::clipped_cell(size_t index, Inset inset) const {
    CellClipper cell = fundamental_region(index, inset);
    Vector2 own = site(index);

    for (const Neighbor& neighbor : _triangulation.neighbors(index)) {
        if (neighbor.index == index) {
            continue;
        }

        Vector2 inward = own - neighbor.position;

        cell.clip(
            inward,
            0.5 * (own + neighbor.position) + inset.bisector * inward.normalized(),
            Bisector{neighbor.index, neighbor.position}
        );
    }

    return cell;
}

// The rectangle every point of the cell must lie in: along a wrapped axis
// the strip closer to the site than to its period images, along a walled
// axis the domain itself, each side moved by the inset its cut takes.
inline CellClipper Diagram::fundamental_region(size_t index, Inset inset) const {
    Vector2 own = site(index);
    Vector2 low;
    Vector2 high;
    std::vector<Cut> low_cuts;
    std::vector<Cut> high_cuts;

    for (int axis = 0; axis < 2; ++axis) {
        double extent = domain().extent(axis);
        Vector2 period = Vector2::Zero();
        period[axis] = extent;
        Vector2 inward = Vector2::Zero();
        inward[axis] = 1.0;

        if (domain().wrapped(axis)) {
            low[axis] = own[axis] - 0.5 * extent + inset.bisector;
            high[axis] = own[axis] + 0.5 * extent - inset.bisector;
            low_cuts.push_back(Bisector{index, own - period});
            high_cuts.push_back(Bisector{index, own + period});
        } else {
            low[axis] = inset.wall;
            high[axis] = extent - inset.wall;
            low_cuts.push_back(Wall{inward});
            high_cuts.push_back(Wall{-inward});
        }
    }

    if (low.x() >= high.x() || low.y() >= high.y()) {
        return CellClipper(std::vector<CellClipper::Edge>{});
    }

    return CellClipper(std::vector<CellClipper::Edge>{
        CellClipper::Edge{Vector2(low.x(), low.y()), low_cuts[1]},
        CellClipper::Edge{Vector2(high.x(), low.y()), high_cuts[0]},
        CellClipper::Edge{Vector2(high.x(), high.y()), high_cuts[1]},
        CellClipper::Edge{Vector2(low.x(), high.y()), low_cuts[0]}
    });
}

} // namespace geometry_art::voronoi::flat

#endif //GEOMETRY_ART_VORONOI_FLAT_CORE_DIAGRAM_HPP_
