#ifndef GEOMETRY_ART_VORONOI_EDGE_SLOTS_HPP_
#define GEOMETRY_ART_VORONOI_EDGE_SLOTS_HPP_

#include <cstddef>
#include <map>
#include <utility>
#include <vector>

namespace geometry_art::voronoi {

// One slot per unique bisector. Each bisector borders two cells and appears
// in both edge lists; anything integrated along it is computed once, in the
// slot, and read from either side.
struct EdgeSlots {
    std::vector<std::vector<size_t>> slots_by_cell;

    // One (cell, position) pair per slot, naming an edge that carries the
    // slot's geometry.
    std::vector<std::pair<size_t, size_t>> representatives;

    [[nodiscard]] size_t count() const { return representatives.size(); }

    template<typename EdgeType>
    [[nodiscard]] static EdgeSlots build(const std::vector<std::vector<EdgeType>>& edges_by_cell);
};

template<typename EdgeType>
EdgeSlots EdgeSlots::build(const std::vector<std::vector<EdgeType>>& edges_by_cell) {
    using EdgeKey = std::pair<size_t, size_t>;

    std::map<EdgeKey, size_t> slot_by_key;
    EdgeSlots slots;
    slots.slots_by_cell.resize(edges_by_cell.size());

    for (size_t cell = 0; cell < edges_by_cell.size(); ++cell) {
        slots.slots_by_cell[cell].reserve(edges_by_cell[cell].size());

        for (size_t position = 0; position < edges_by_cell[cell].size(); ++position) {
            size_t neighbor = edges_by_cell[cell][position].neighbor_index;
            EdgeKey key = cell < neighbor ? EdgeKey{cell, neighbor} : EdgeKey{neighbor, cell};
            auto [iterator, inserted] = slot_by_key.try_emplace(key, slots.representatives.size());

            if (inserted) {
                slots.representatives.emplace_back(cell, position);
            }

            slots.slots_by_cell[cell].push_back(iterator->second);
        }
    }

    return slots;
}

} // namespace geometry_art::voronoi

#endif //GEOMETRY_ART_VORONOI_EDGE_SLOTS_HPP_
