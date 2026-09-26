#ifndef GEOMETRY_ART_IO_SNAPSHOT_JSON_READER_HPP_
#define GEOMETRY_ART_IO_SNAPSHOT_JSON_READER_HPP_

#include "snapshot.hpp"
#include "../../types.hpp"
#include <simdjson.h>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <string_view>

namespace geometry_art::io::snapshot {

// Reads a snapshot back from the JSON the writer emits, so anything
// downstream of the solver can start from a saved tessellation.
class JsonReader {
 public:
    [[nodiscard]] Snapshot read(const std::string& json) const;
    [[nodiscard]] Snapshot read_file(const std::filesystem::path& path) const;

 private:
    [[nodiscard]] static Snapshot::Cell read_cell(simdjson::ondemand::object cell);
    [[nodiscard]] static Vector3 read_point(simdjson::ondemand::array point);
};

inline Snapshot JsonReader::read(const std::string& json) const {
    simdjson::padded_string padded(json);
    simdjson::ondemand::parser parser;

    try {
        simdjson::ondemand::document document = parser.iterate(padded);
        Snapshot snapshot;
        snapshot.geometry = std::string(std::string_view(document["geometry"]));
        snapshot.width = double(document["width"]);
        snapshot.height = double(document["height"]);
        double scale = 0.0;

        if (document["scale"].get_double().get(scale) == simdjson::SUCCESS) {
            snapshot.scale = scale;
        }

        snapshot.total_mass = double(document["totalMass"]);

        for (simdjson::ondemand::object cell : document["cells"]) {
            snapshot.cells.push_back(read_cell(cell));
        }

        return snapshot;
    } catch (const simdjson::simdjson_error& error) {
        throw std::runtime_error(std::string("not a snapshot: ") + error.what());
    }
}

inline Snapshot JsonReader::read_file(const std::filesystem::path& path) const {
    simdjson::padded_string padded;

    if (simdjson::padded_string::load(path.string()).get(padded) != simdjson::SUCCESS) {
        throw std::runtime_error("cannot read " + path.string());
    }

    return read(std::string(padded.data(), padded.size()));
}

inline Snapshot::Cell JsonReader::read_cell(simdjson::ondemand::object cell) {
    Snapshot::Cell entry;
    entry.site_index = uint64_t(cell["site"]);
    entry.site = read_point(cell["position"]);
    entry.mass = double(cell["mass"]);
    entry.area = double(cell["area"]);

    for (uint64_t neighbor : cell["neighbors"]) {
        entry.neighbors.push_back(neighbor);
    }

    for (simdjson::ondemand::array point : cell["boundary"]) {
        entry.boundary.push_back(read_point(point));
    }

    return entry;
}

inline Vector3 JsonReader::read_point(simdjson::ondemand::array point) {
    Vector3 result = Vector3::Zero();
    int axis = 0;

    for (double coordinate : point) {
        if (axis < 3) {
            result[axis] = coordinate;
        }

        ++axis;
    }

    return result;
}

} // namespace geometry_art::io::snapshot

#endif //GEOMETRY_ART_IO_SNAPSHOT_JSON_READER_HPP_
