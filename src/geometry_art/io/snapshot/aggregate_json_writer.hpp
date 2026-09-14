#ifndef GEOMETRY_ART_IO_SNAPSHOT_AGGREGATE_JSON_WRITER_HPP_
#define GEOMETRY_ART_IO_SNAPSHOT_AGGREGATE_JSON_WRITER_HPP_

#include "aggregate_snapshot.hpp"
#include "../../types.hpp"
#include <iomanip>
#include <ostream>
#include <sstream>
#include <string>

namespace geometry_art::io::snapshot {

// Writes an aggregate snapshot as JSON with fixed digits, for the same reason
// as JsonWriter: a golden file that reformats on a different platform is not
// a golden file.
class AggregateJsonWriter {
 public:
    explicit AggregateJsonWriter(int precision = 9);

    void write(const AggregateSnapshot& snapshot, std::ostream& stream) const;
    [[nodiscard]] std::string to_string(const AggregateSnapshot& snapshot) const;

 private:
    int _precision;

    void write_particle(const dla::Particle& particle, std::ostream& stream) const;
    void write_point(const Vector3& point, std::ostream& stream) const;
    [[nodiscard]] std::string number(double value) const;
};

inline AggregateJsonWriter::AggregateJsonWriter(int precision) :
    _precision(precision) {
}

inline std::string AggregateJsonWriter::to_string(const AggregateSnapshot& snapshot) const {
    std::ostringstream stream;
    write(snapshot, stream);
    return stream.str();
}

inline void AggregateJsonWriter::write(const AggregateSnapshot& snapshot, std::ostream& stream) const {
    stream << "{\n";
    stream << "  \"particleRadius\": " << number(snapshot.particle_radius) << ",\n";
    stream << "  \"reach\": " << number(snapshot.reach) << ",\n";
    stream << "  \"particles\": [\n";

    for (size_t index = 0; index < snapshot.particles.size(); ++index) {
        write_particle(snapshot.particles[index], stream);
        stream << (index + 1 < snapshot.particles.size() ? ",\n" : "\n");
    }

    stream << "  ]\n";
    stream << "}\n";
}

inline void AggregateJsonWriter::write_particle(const dla::Particle& particle, std::ostream& stream) const {
    stream << "    {\"center\": ";
    write_point(particle.center, stream);
    stream << ", \"parent\": ";

    if (particle.parent.has_value()) {
        stream << *particle.parent;
    } else {
        stream << "null";
    }

    stream << "}";
}

inline void AggregateJsonWriter::write_point(const Vector3& point, std::ostream& stream) const {
    stream << "[" << number(point.x()) << ", " << number(point.y()) << ", " << number(point.z()) << "]";
}

inline std::string AggregateJsonWriter::number(double value) const {
    std::ostringstream stream;
    stream << std::setprecision(_precision) << (value == 0.0 ? 0.0 : value);
    return stream.str();
}

} // namespace geometry_art::io::snapshot

#endif //GEOMETRY_ART_IO_SNAPSHOT_AGGREGATE_JSON_WRITER_HPP_
