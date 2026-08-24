#ifndef GLOBEART_SRC_GLOBE_IO_SNAPSHOT_JSON_WRITER_HPP_
#define GLOBEART_SRC_GLOBE_IO_SNAPSHOT_JSON_WRITER_HPP_

#include "snapshot.hpp"
#include <iomanip>
#include <ostream>
#include <sstream>
#include <string>

namespace globe::io::snapshot {

// Writes a snapshot as JSON for the viewer and for golden files, so the
// digits are fixed rather than left to the stream's defaults: a golden file
// that reformats on a different platform is not a golden file.
class JsonWriter {
 public:
    explicit JsonWriter(int precision = 9);

    void write(const Snapshot& snapshot, std::ostream& stream) const;
    [[nodiscard]] std::string to_string(const Snapshot& snapshot) const;

 private:
    int _precision;

    void write_cell(const Snapshot::Cell& cell, std::ostream& stream) const;
    void write_point(const Vector3& point, std::ostream& stream) const;
    [[nodiscard]] std::string number(double value) const;
};

inline JsonWriter::JsonWriter(int precision) :
    _precision(precision) {
}

inline std::string JsonWriter::to_string(const Snapshot& snapshot) const {
    std::ostringstream stream;
    write(snapshot, stream);
    return stream.str();
}

inline void JsonWriter::write(const Snapshot& snapshot, std::ostream& stream) const {
    stream << "{\n";
    stream << "  \"geometry\": \"" << snapshot.geometry << "\",\n";
    stream << "  \"width\": " << number(snapshot.width) << ",\n";
    stream << "  \"height\": " << number(snapshot.height) << ",\n";
    stream << "  \"totalMass\": " << number(snapshot.total_mass) << ",\n";
    stream << "  \"targetMass\": " << number(snapshot.target_mass()) << ",\n";
    stream << "  \"relativeRmsCapacityError\": " << number(snapshot.relative_rms_capacity_error()) << ",\n";
    stream << "  \"cells\": [\n";

    for (size_t index = 0; index < snapshot.cells.size(); ++index) {
        write_cell(snapshot.cells[index], stream);
        stream << (index + 1 < snapshot.cells.size() ? ",\n" : "\n");
    }

    stream << "  ]\n";
    stream << "}\n";
}

inline void JsonWriter::write_cell(const Snapshot::Cell& cell, std::ostream& stream) const {
    stream << "    {\"site\": " << cell.site_index << ", \"position\": ";
    write_point(cell.site, stream);
    stream << ", \"mass\": " << number(cell.mass);
    stream << ", \"area\": " << number(cell.area);
    stream << ", \"neighbors\": [";

    for (size_t index = 0; index < cell.neighbors.size(); ++index) {
        stream << (index == 0 ? "" : ", ") << cell.neighbors[index];
    }

    stream << "], \"boundary\": [";

    for (size_t index = 0; index < cell.boundary.size(); ++index) {
        stream << (index == 0 ? "" : ", ");
        write_point(cell.boundary[index], stream);
    }

    stream << "]}";
}

inline void JsonWriter::write_point(const Vector3& point, std::ostream& stream) const {
    stream << "[" << number(point.x()) << ", " << number(point.y()) << ", " << number(point.z()) << "]";
}

// Negative zero and a stray exponent format are the two ways an otherwise
// identical value produces a differing golden file.
inline std::string JsonWriter::number(double value) const {
    std::ostringstream stream;
    stream << std::setprecision(_precision) << (value == 0.0 ? 0.0 : value);
    return stream.str();
}

} // namespace globe::io::snapshot

#endif //GLOBEART_SRC_GLOBE_IO_SNAPSHOT_JSON_WRITER_HPP_
