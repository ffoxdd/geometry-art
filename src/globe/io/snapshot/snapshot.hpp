#ifndef GLOBEART_SRC_GLOBE_IO_SNAPSHOT_SNAPSHOT_HPP_
#define GLOBEART_SRC_GLOBE_IO_SNAPSHOT_SNAPSHOT_HPP_

#include "../../types.hpp"
#include "../../fields/spherical/field.hpp"
#include "../../std_ext/parallel_for.hpp"
#include "../../voronoi/spherical/core/diagram.hpp"
#include <cstddef>
#include <string>
#include <vector>

namespace globe::io::snapshot {

using voronoi::spherical::CellEdgeInfo;
using voronoi::spherical::Diagram;

// A tessellation as plain data: no triangulation, no field, no renderer.
// Everything downstream of the solver reads this -- the viewer, the vector
// output, the golden files -- so what the pipeline produces can be compared
// and drawn without linking any of it.
struct Snapshot {
    struct Cell {
        size_t site_index = 0;
        Vector3 site;
        std::vector<Vector3> boundary;
        std::vector<size_t> neighbors;
        double mass = 0.0;
        double area = 0.0;
    };

    std::string geometry;

    // The rectangle of periods, for flat geometries; unused on the sphere.
    double width = 0.0;
    double height = 0.0;

    std::vector<Cell> cells;
    double total_mass = 0.0;

    [[nodiscard]] double target_mass() const;
    [[nodiscard]] double relative_rms_capacity_error() const;
};

template<Diagram DiagramType, fields::spherical::Field FieldType>
[[nodiscard]] Snapshot capture(const DiagramType& diagram, const FieldType& field);

inline double Snapshot::target_mass() const {
    if (cells.empty()) {
        return 0.0;
    }

    return total_mass / static_cast<double>(cells.size());
}

inline double Snapshot::relative_rms_capacity_error() const {
    double target = target_mass();

    if (cells.empty() || target == 0.0) {
        return 0.0;
    }

    double sum = 0.0;

    for (const Cell& cell : cells) {
        double error = cell.mass - target;
        sum += error * error;
    }

    return std::sqrt(sum / static_cast<double>(cells.size())) / target;
}

template<Diagram DiagramType, fields::spherical::Field FieldType>
Snapshot capture(const DiagramType& diagram, const FieldType& field) {
    Snapshot snapshot;
    snapshot.geometry = "sphere";
    snapshot.total_mass = field.total_mass();
    snapshot.cells.resize(diagram.size());

    std_ext::parallel_for(diagram.size(), [&](size_t index) {
        Polygon cell = diagram.cell(index);
        Snapshot::Cell& entry = snapshot.cells[index];
        entry.site_index = index;
        entry.site = to_vector3(diagram.site(index));
        entry.mass = field.integrals(cell).mass;
        entry.area = cell.area();

        for (const Arc& arc : cell.arcs()) {
            entry.boundary.push_back(Vector3(arc.source()));
        }

        for (const CellEdgeInfo& edge : diagram.cell_edges(index)) {
            entry.neighbors.push_back(edge.neighbor_index);
        }
    });

    return snapshot;
}

} // namespace globe::io::snapshot

#endif //GLOBEART_SRC_GLOBE_IO_SNAPSHOT_SNAPSHOT_HPP_
