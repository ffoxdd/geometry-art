#ifndef GEOMETRY_ART_IO_SNAPSHOT_AGGREGATE_SVG_WRITER_HPP_
#define GEOMETRY_ART_IO_SNAPSHOT_AGGREGATE_SVG_WRITER_HPP_

#include "aggregate_snapshot.hpp"
#include "../../types.hpp"
#include <algorithm>
#include <iomanip>
#include <numeric>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

namespace geometry_art::io::snapshot {

// An orthographic view of an aggregate as SVG: one circle per particle,
// painted back to front and shaded by depth, scaled so the aggregate's reach
// fills the frame.
class AggregateSvgWriter {
 public:
    explicit AggregateSvgWriter(double size = 512.0);

    void write(const AggregateSnapshot& snapshot, std::ostream& stream) const;
    [[nodiscard]] std::string to_string(const AggregateSnapshot& snapshot) const;

 private:
    double _size;

    void write_particle(const AggregateSnapshot& snapshot, const dla::Particle& particle, std::ostream& stream) const;
    [[nodiscard]] static std::vector<size_t> back_to_front(const AggregateSnapshot& snapshot);
    [[nodiscard]] std::string shade(const AggregateSnapshot& snapshot, const dla::Particle& particle) const;
    [[nodiscard]] double scale(const AggregateSnapshot& snapshot) const;
    [[nodiscard]] std::string number(double value) const;
};

inline AggregateSvgWriter::AggregateSvgWriter(double size) :
    _size(size) {
}

inline std::string AggregateSvgWriter::to_string(const AggregateSnapshot& snapshot) const {
    std::ostringstream stream;
    write(snapshot, stream);
    return stream.str();
}

inline void AggregateSvgWriter::write(const AggregateSnapshot& snapshot, std::ostream& stream) const {
    stream << "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 "
        << number(_size) << " " << number(_size) << "\">\n";
    stream << "  <rect width=\"" << number(_size) << "\" height=\"" << number(_size)
        << "\" fill=\"#f7f5f2\"/>\n";

    for (size_t index : back_to_front(snapshot)) {
        write_particle(snapshot, snapshot.particles[index], stream);
    }

    stream << "</svg>\n";
}

inline void AggregateSvgWriter::write_particle(
    const AggregateSnapshot& snapshot,
    const dla::Particle& particle,
    std::ostream& stream
) const {
    double centre = _size / 2.0;

    stream << "  <circle cx=\"" << number(centre + scale(snapshot) * particle.center.x())
        << "\" cy=\"" << number(centre - scale(snapshot) * particle.center.y())
        << "\" r=\"" << number(scale(snapshot) * snapshot.particle_radius)
        << "\" fill=\"" << shade(snapshot, particle) << "\"/>\n";
}

inline std::vector<size_t> AggregateSvgWriter::back_to_front(const AggregateSnapshot& snapshot) {
    std::vector<size_t> order(snapshot.particles.size());
    std::iota(order.begin(), order.end(), 0);

    std::stable_sort(order.begin(), order.end(), [&snapshot](size_t a, size_t b) {
        return snapshot.particles[a].center.z() < snapshot.particles[b].center.z();
    });

    return order;
}

// Depth as brightness: the far side fades toward the background, the near
// side stays dark, so the flat projection still reads as a volume.
inline std::string AggregateSvgWriter::shade(
    const AggregateSnapshot& snapshot,
    const dla::Particle& particle
) const {
    double depth = (particle.center.z() + snapshot.reach) / (2.0 * snapshot.reach);
    int brightness = static_cast<int>(200.0 - 170.0 * std::clamp(depth, 0.0, 1.0));

    std::ostringstream color;
    color << "rgb(" << brightness << "," << brightness << "," << brightness << ")";
    return color.str();
}

inline double AggregateSvgWriter::scale(const AggregateSnapshot& snapshot) const {
    return (_size * 0.48) / snapshot.reach;
}

inline std::string AggregateSvgWriter::number(double value) const {
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(2) << (value == 0.0 ? 0.0 : value);
    return stream.str();
}

} // namespace geometry_art::io::snapshot

#endif //GEOMETRY_ART_IO_SNAPSHOT_AGGREGATE_SVG_WRITER_HPP_
