#ifndef GEOMETRY_ART_GEOMETRY_PLANAR_DOMAIN_DELAUNAY_HPP_
#define GEOMETRY_ART_GEOMETRY_PLANAR_DOMAIN_DELAUNAY_HPP_

#include "domain.hpp"
#include "../../types.hpp"
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/assertions.h>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

namespace geometry_art::geometry::planar {

// A Delaunay neighbour of a site, placed in that site's own chart: a
// neighbour across a wrapped seam reads as the period image that actually
// borders it, so no caller handles an offset.
struct Neighbor {
    size_t index;
    Vector2 position;
};

// A Delaunay triangulation of sites on a flat domain: periodic along each
// wrapped axis, bounded along each walled one. A site's Delaunay neighbours
// are exactly the sites whose bisectors can bound its Voronoi cell, which is
// all the diagram reads from here.
//
// Wrapping is computed literally, by tiling the sites along the wrapped
// axes and reading each site's neighbours off its own copy. Interior sites
// need no copies at all -- a copy can only be a Delaunay neighbour of the
// domain's own sites if it lies within the largest circumradius of the
// domain, which is a few cell widths rather than a period -- so images are
// laid in a band around the seam whose width is measured rather than
// assumed. A walled axis lays no images: its hull is the domain's edge, and
// the infinite faces there are the triangulation's honest answer.
class DomainDelaunay {
 public:
    DomainDelaunay(Domain domain, std::vector<Vector2> sites);

    DomainDelaunay(const DomainDelaunay&) = delete;
    DomainDelaunay& operator=(const DomainDelaunay&) = delete;
    DomainDelaunay(DomainDelaunay&&) = default;
    DomainDelaunay& operator=(DomainDelaunay&&) = default;

    [[nodiscard]] size_t size() const { return _handles.size(); }
    [[nodiscard]] const Domain& domain() const { return _domain; }

    [[nodiscard]] Vector2 site(size_t index) const;
    [[nodiscard]] std::vector<Neighbor> neighbors(size_t index) const;

 private:
    // A periodic image carries the index of the site it copies, so a wrapped
    // neighbour reports the same index its canonical placement would, and a
    // flag for the one copy that lies in the domain and owns the cell.
    struct VertexInfo {
        size_t site_index;
        bool central;
    };

    using Kernel = ::CGAL::Exact_predicates_inexact_constructions_kernel;
    using VertexBase = ::CGAL::Triangulation_vertex_base_with_info_2<VertexInfo, Kernel>;
    using DataStructure = ::CGAL::Triangulation_data_structure_2<VertexBase>;
    using Triangulation = ::CGAL::Delaunay_triangulation_2<Kernel, DataStructure>;
    using VertexHandle = Triangulation::Vertex_handle;
    using FaceHandle = Triangulation::Face_handle;

    // How many cell widths of images to lay around the domain before the
    // covering is measured. Wide enough that one pass almost always suffices,
    // and cheap enough that overshooting costs little.
    static constexpr double BAND_CELL_LAYERS = 4.0;

    Domain _domain;
    std::unique_ptr<Triangulation> _triangulation;
    std::vector<VertexHandle> _handles;

    void build(const std::vector<Vector2>& sites, double band);
    [[nodiscard]] double initial_band(size_t count) const;
    [[nodiscard]] double shortest_period() const;
    [[nodiscard]] double distance_to_domain(const Vector2& point) const;
    [[nodiscard]] double maximum_circumradius() const;
    [[nodiscard]] double circumradius(FaceHandle face) const;
    [[nodiscard]] Vector2 position(VertexHandle vertex) const;
};

// The band is laid at a guess, measured against the circumradii it produced,
// and widened until it covers them -- terminating at the full three-tile
// replication, which is the most a covering can be and which suffices exactly
// when no circumradius reaches half a period.
inline DomainDelaunay::DomainDelaunay(Domain domain, std::vector<Vector2> sites) :
    _domain(domain),
    _triangulation(std::make_unique<Triangulation>()) {
    for (Vector2& site : sites) {
        site = _domain.canonical(site);
    }

    if (!_domain.has_wrapped_axis()) {
        build(sites, 0.0);
        return;
    }

    double full_replication = std::max(_domain.width, _domain.height);
    double band = initial_band(sites.size());

    while (true) {
        build(sites, band);
        double radius = maximum_circumradius();

        if (radius < band || band >= full_replication) {
            CGAL_assertion(radius < 0.5 * shortest_period());
            static_cast<void>(radius);
            return;
        }

        band = std::min(2.0 * band, full_replication);
    }
}

inline Vector2 DomainDelaunay::site(size_t index) const {
    return position(_handles[index]);
}

// A lone site has no fan to walk: the triangulation is a single vertex.
inline std::vector<Neighbor> DomainDelaunay::neighbors(size_t index) const {
    std::vector<Neighbor> result;

    if (_triangulation->dimension() < 1) {
        return result;
    }

    VertexHandle vertex = _handles[index];
    auto circulator = _triangulation->incident_vertices(vertex);
    auto start = circulator;

    do {
        VertexHandle neighbor = circulator;

        if (!_triangulation->is_infinite(neighbor)) {
            result.push_back(Neighbor{neighbor->info().site_index, position(neighbor)});
        }

        ++circulator;
    } while (circulator != start);

    return result;
}

inline double DomainDelaunay::initial_band(size_t count) const {
    if (count == 0) {
        return std::max(_domain.width, _domain.height);
    }

    return BAND_CELL_LAYERS * std::sqrt(_domain.area() / static_cast<double>(count));
}

inline double DomainDelaunay::shortest_period() const {
    double shortest = std::numeric_limits<double>::infinity();

    for (int axis = 0; axis < 2; ++axis) {
        if (_domain.wrapped(axis)) {
            shortest = std::min(shortest, _domain.extent(axis));
        }
    }

    return shortest;
}

inline void DomainDelaunay::build(const std::vector<Vector2>& sites, double band) {
    _triangulation = std::make_unique<Triangulation>();
    _handles.assign(sites.size(), VertexHandle());

    int column_reach = _domain.wrapped(0) ? 1 : 0;
    int row_reach = _domain.wrapped(1) ? 1 : 0;

    std::vector<std::pair<Triangulation::Point, VertexInfo>> points;
    points.reserve(sites.size());

    for (size_t index = 0; index < sites.size(); ++index) {
        for (int row = -row_reach; row <= row_reach; ++row) {
            for (int column = -column_reach; column <= column_reach; ++column) {
                bool central = row == 0 && column == 0;
                Vector2 copy = sites[index] + Vector2(column * _domain.width, row * _domain.height);

                if (!central && distance_to_domain(copy) > band) {
                    continue;
                }

                points.emplace_back(
                    Triangulation::Point(copy.x(), copy.y()),
                    VertexInfo{index, central}
                );
            }
        }
    }

    // The range insertion spatially sorts the points, so each location walk
    // starts beside the last instead of crossing the triangulation.
    _triangulation->insert(points.begin(), points.end());

    for (auto vertex = _triangulation->finite_vertices_begin();
        vertex != _triangulation->finite_vertices_end(); ++vertex) {
        if (vertex->info().central) {
            _handles[vertex->info().site_index] = vertex;
        }
    }
}

inline double DomainDelaunay::distance_to_domain(const Vector2& point) const {
    double across = std::max({0.0, -point.x(), point.x() - _domain.width});
    double along = std::max({0.0, -point.y(), point.y() - _domain.height});

    return std::sqrt(across * across + along * along);
}

// A site whose own fan touches the edge of the covered region has neighbours
// the band left out, which reads as an infinite incident face. Where the
// domain has a wall, the hull along that wall is genuine: the tiling never
// extends past it, so an infinite face there says nothing about the band.
inline double DomainDelaunay::maximum_circumradius() const {
    double maximum = 0.0;

    for (const VertexHandle& vertex : _handles) {
        if (vertex == VertexHandle()) {
            continue;
        }

        auto circulator = _triangulation->incident_faces(vertex);
        auto start = circulator;

        do {
            if (_triangulation->is_infinite(circulator)) {
                if (!_domain.has_walled_axis()) {
                    return std::numeric_limits<double>::infinity();
                }
            } else {
                maximum = std::max(maximum, circumradius(circulator));
            }

            ++circulator;
        } while (circulator != start);
    }

    return maximum;
}

inline double DomainDelaunay::circumradius(FaceHandle face) const {
    auto center = _triangulation->dual(face);
    auto corner = _triangulation->point(face->vertex(0));

    return std::sqrt(::CGAL::to_double(::CGAL::squared_distance(center, corner)));
}

inline Vector2 DomainDelaunay::position(VertexHandle vertex) const {
    auto point = _triangulation->point(vertex);
    return Vector2(point.x(), point.y());
}

} // namespace geometry_art::geometry::planar

#endif //GEOMETRY_ART_GEOMETRY_PLANAR_DOMAIN_DELAUNAY_HPP_
