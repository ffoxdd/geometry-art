#ifndef GEOMETRY_ART_IO_SNAPSHOT_SVG_WRITER_HPP_
#define GEOMETRY_ART_IO_SNAPSHOT_SVG_WRITER_HPP_

#include "snapshot.hpp"
#include "../../types.hpp"
#include <cmath>
#include <iomanip>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

namespace geometry_art::io::snapshot {

// An orthographic view of a snapshot as SVG. Being text, it diffs in review
// and needs no browser or pixel comparison to serve as a regression check on
// what the tessellation looks like, not just on what it measures.
class SvgWriter {
 public:
    explicit SvgWriter(double size = 512.0, int samples_per_edge = 8);

    void write(const Snapshot& snapshot, std::ostream& stream) const;
    [[nodiscard]] std::string to_string(const Snapshot& snapshot) const;

 private:
    double _size;
    int _samples_per_edge;

    void write_cell(const Snapshot::Cell& cell, std::ostream& stream) const;
    [[nodiscard]] std::vector<Vector3> boundary_samples(const Snapshot::Cell& cell) const;
    [[nodiscard]] std::string path_of(const std::vector<Vector3>& samples) const;
    [[nodiscard]] std::string projected(const Vector3& point) const;
    [[nodiscard]] std::string number(double value) const;
};

inline SvgWriter::SvgWriter(double size, int samples_per_edge) :
    _size(size),
    _samples_per_edge(samples_per_edge) {
}

inline std::string SvgWriter::to_string(const Snapshot& snapshot) const {
    std::ostringstream stream;
    write(snapshot, stream);
    return stream.str();
}

inline void SvgWriter::write(const Snapshot& snapshot, std::ostream& stream) const {
    stream << "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 "
        << number(_size) << " " << number(_size) << "\">\n";
    stream << "  <circle cx=\"" << number(_size / 2.0) << "\" cy=\"" << number(_size / 2.0)
        << "\" r=\"" << number(_size * 0.48) << "\" fill=\"#f7f5f2\"/>\n";

    for (const Snapshot::Cell& cell : snapshot.cells) {
        write_cell(cell, stream);
    }

    stream << "</svg>\n";
}

// Only the hemisphere facing the viewer is drawn, so a cell that straddles
// the silhouette contributes the part of its boundary that is in front.
inline void SvgWriter::write_cell(const Snapshot::Cell& cell, std::ostream& stream) const {
    if (cell.site.z() <= 0.0) {
        return;
    }

    std::vector<Vector3> samples = boundary_samples(cell);
    std::string path = path_of(samples);

    if (path.empty()) {
        return;
    }

    stream << "  <path d=\"" << path << "\" fill=\"none\" stroke=\"#1c1c1c\" stroke-width=\"1.2\"/>\n";
}

inline std::vector<Vector3> SvgWriter::boundary_samples(const Snapshot::Cell& cell) const {
    std::vector<Vector3> samples;

    for (size_t index = 0; index < cell.boundary.size(); ++index) {
        const Vector3& from = cell.boundary[index];
        const Vector3& to = cell.boundary[(index + 1) % cell.boundary.size()];

        for (int step = 0; step < _samples_per_edge; ++step) {
            double fraction = static_cast<double>(step) / _samples_per_edge;
            Vector3 blended = from + fraction * (to - from);
            double norm = blended.norm();
            samples.push_back(norm > 0.0 ? Vector3(blended / norm) : from);
        }
    }

    return samples;
}

inline std::string SvgWriter::path_of(const std::vector<Vector3>& samples) const {
    std::ostringstream path;
    bool started = false;

    for (const Vector3& sample : samples) {
        if (sample.z() < 0.0) {
            started = false;
            continue;
        }

        path << (started ? " L " : (path.tellp() > 0 ? " M " : "M ")) << projected(sample);
        started = true;
    }

    return path.str();
}

inline std::string SvgWriter::projected(const Vector3& point) const {
    double radius = _size * 0.48;
    double centre = _size / 2.0;
    return number(centre + radius * point.x()) + " " + number(centre - radius * point.y());
}

inline std::string SvgWriter::number(double value) const {
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(2) << (value == 0.0 ? 0.0 : value);
    return stream.str();
}

} // namespace geometry_art::io::snapshot

#endif //GEOMETRY_ART_IO_SNAPSHOT_SVG_WRITER_HPP_
