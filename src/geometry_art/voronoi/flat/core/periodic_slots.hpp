#ifndef GEOMETRY_ART_VORONOI_FLAT_CORE_PERIODIC_SLOTS_HPP_
#define GEOMETRY_ART_VORONOI_FLAT_CORE_PERIODIC_SLOTS_HPP_

#include "diagram.hpp"
#include <cmath>
#include <cstddef>
#include <map>
#include <tuple>
#include <utility>
#include <vector>

namespace geometry_art::voronoi::flat {

// One slot per physical bisector. On a wrapped domain the unordered site pair is not
// enough to name one: the same two sites can share two bisectors, one
// across each seam, and a cell can border itself. The key adds the period
// offset between the two charts, oriented from the smaller index so both
// sides produce the same key; a self-edge's two appearances get opposite
// offsets and therefore separate slots, which keeps their cancelling sweep
// contributions independently accounted.
struct PeriodicSlots {
    std::vector<std::vector<size_t>> slots_by_cell;
    std::vector<std::pair<size_t, size_t>> representatives;

    [[nodiscard]] size_t count() const { return representatives.size(); }

    [[nodiscard]] static PeriodicSlots build(
        const Diagram& diagram,
        const std::vector<std::vector<CellEdgeInfo>>& edges_by_cell
    );

 private:
    using Key = std::tuple<size_t, size_t, int, int>;

    [[nodiscard]] static Key key_of(const Diagram& diagram, size_t cell, const CellEdgeInfo& edge);
};

inline PeriodicSlots PeriodicSlots::build(
    const Diagram& diagram,
    const std::vector<std::vector<CellEdgeInfo>>& edges_by_cell
) {
    std::map<Key, size_t> slot_by_key;
    PeriodicSlots slots;
    slots.slots_by_cell.resize(edges_by_cell.size());

    for (size_t cell = 0; cell < edges_by_cell.size(); ++cell) {
        slots.slots_by_cell[cell].reserve(edges_by_cell[cell].size());

        for (size_t position = 0; position < edges_by_cell[cell].size(); ++position) {
            auto [iterator, inserted] = slot_by_key.try_emplace(
                key_of(diagram, cell, edges_by_cell[cell][position]),
                slots.representatives.size()
            );

            if (inserted) {
                slots.representatives.emplace_back(cell, position);
            }

            slots.slots_by_cell[cell].push_back(iterator->second);
        }
    }

    return slots;
}

inline PeriodicSlots::Key PeriodicSlots::key_of(const Diagram& diagram, size_t cell, const CellEdgeInfo& edge) {
    Vector2 offset = edge.neighbor_position - diagram.site(edge.neighbor_index);
    int dx = static_cast<int>(std::lround(offset.x() / diagram.width()));
    int dy = static_cast<int>(std::lround(offset.y() / diagram.height()));

    if (cell <= edge.neighbor_index) {
        return Key{cell, edge.neighbor_index, dx, dy};
    }

    return Key{edge.neighbor_index, cell, -dx, -dy};
}

} // namespace geometry_art::voronoi::flat

#endif //GEOMETRY_ART_VORONOI_FLAT_CORE_PERIODIC_SLOTS_HPP_
