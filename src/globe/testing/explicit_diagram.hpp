#ifndef GLOBEART_SRC_GLOBE_TESTING_EXPLICIT_DIAGRAM_HPP_
#define GLOBEART_SRC_GLOBE_TESTING_EXPLICIT_DIAGRAM_HPP_

#include "../voronoi/spherical/core/diagram.hpp"
#include <cstddef>
#include <utility>
#include <vector>

namespace globe::testing {

using voronoi::spherical::CellEdgeInfo;
using voronoi::spherical::Polygon;

// A tessellation given cell by cell instead of derived from sites. It
// satisfies the same Diagram concept the solver reads, so a scenario with two
// or three cells of known geometry can drive the integration and derivative
// machinery without a triangulation deciding what the cells are.
class ExplicitDiagram {
 public:
    void add(cgal::Point3 site, Polygon cell, std::vector<CellEdgeInfo> edges);

    [[nodiscard]] size_t size() const { return _sites.size(); }
    [[nodiscard]] cgal::Point3 site(size_t index) const { return _sites[index]; }
    [[nodiscard]] Polygon cell(size_t index) const { return _cells[index]; }
    [[nodiscard]] std::vector<CellEdgeInfo> cell_edges(size_t index) const { return _edges[index]; }

    // The two hemispheres of a pair of antipodal sites: the smallest
    // tessellation whose cells, masses and shared bisector are all known in
    // closed form.
    [[nodiscard]] static ExplicitDiagram hemispheres();

 private:
    std::vector<cgal::Point3> _sites;
    std::vector<Polygon> _cells;
    std::vector<std::vector<CellEdgeInfo>> _edges;
};

inline void ExplicitDiagram::add(cgal::Point3 site, Polygon cell, std::vector<CellEdgeInfo> edges) {
    _sites.push_back(site);
    _cells.push_back(std::move(cell));
    _edges.push_back(std::move(edges));
}

inline ExplicitDiagram ExplicitDiagram::hemispheres() {
    VectorS2 x(1, 0, 0);
    VectorS2 y(0, 1, 0);
    VectorS2 negative_x(-1, 0, 0);
    VectorS2 negative_y(0, -1, 0);
    VectorS2 north(0, 0, 1);

    std::vector<Arc> northern{
        Arc(x, y, north),
        Arc(y, negative_x, north),
        Arc(negative_x, negative_y, north),
        Arc(negative_y, x, north)
    };

    std::vector<Arc> southern{
        Arc(x, negative_y, -north),
        Arc(negative_y, negative_x, -north),
        Arc(negative_x, y, -north),
        Arc(y, x, -north)
    };

    ExplicitDiagram diagram;
    std::vector<CellEdgeInfo> northern_edges;
    std::vector<CellEdgeInfo> southern_edges;

    // Two cells leave no third site at any endpoint, so the opposite
    // indices just name the only neighbour there is.
    for (const Arc& arc : northern) {
        northern_edges.push_back(CellEdgeInfo{1, 1, 1, arc});
    }

    for (const Arc& arc : southern) {
        southern_edges.push_back(CellEdgeInfo{0, 0, 0, arc});
    }

    diagram.add(cgal::to_point(VectorS2(0, 0, 1)), Polygon(northern), std::move(northern_edges));
    diagram.add(cgal::to_point(VectorS2(0, 0, -1)), Polygon(southern), std::move(southern_edges));

    return diagram;
}

} // namespace globe::testing

#endif //GLOBEART_SRC_GLOBE_TESTING_EXPLICIT_DIAGRAM_HPP_
