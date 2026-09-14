#ifndef GEOMETRY_ART_IO_SNAPSHOT_FLAT_CAPTURE_HPP_
#define GEOMETRY_ART_IO_SNAPSHOT_FLAT_CAPTURE_HPP_

#include "snapshot.hpp"
#include "../../fields/flat/field.hpp"
#include "../../std_ext/parallel_for.hpp"
#include "../../voronoi/flat/core/diagram.hpp"
#include <cstddef>
#include <string>

namespace geometry_art::io::snapshot {

// The flat tessellation as the same plain data the sphere produces. Cells
// are reported in their own charts, so on a wrapped axis a cell may
// protrude past the rectangle; the renderer wraps or crops as it pleases.
// The frame size travels in the width and height fields.
template<fields::flat::Field FieldType>
[[nodiscard]] Snapshot capture_flat(
    const voronoi::flat::Diagram& diagram,
    const FieldType& field,
    std::string geometry = "torus"
) {
    Snapshot snapshot;
    snapshot.geometry = std::move(geometry);
    snapshot.width = diagram.width();
    snapshot.height = diagram.height();
    snapshot.total_mass = field.total_mass();
    snapshot.cells.resize(diagram.size());

    std_ext::parallel_for(diagram.size(), [&](size_t index) {
        auto cell = diagram.cell(index);
        Snapshot::Cell& entry = snapshot.cells[index];
        entry.site_index = index;
        entry.site = diagram.site_vector(index);
        entry.mass = field.integrals(cell).mass;
        entry.area = cell.area();

        for (const Vector2& vertex : cell.vertices()) {
            entry.boundary.emplace_back(vertex.x(), vertex.y(), 0.0);
        }

        for (const auto& edge : diagram.cell_edges(index)) {
            entry.neighbors.push_back(edge.neighbor_index);
        }
    });

    return snapshot;
}

} // namespace geometry_art::io::snapshot

#endif //GEOMETRY_ART_IO_SNAPSHOT_FLAT_CAPTURE_HPP_
