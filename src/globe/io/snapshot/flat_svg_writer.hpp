#ifndef GLOBEART_SRC_GLOBE_IO_SNAPSHOT_FLAT_SVG_WRITER_HPP_
#define GLOBEART_SRC_GLOBE_IO_SNAPSHOT_FLAT_SVG_WRITER_HPP_

#include "snapshot.hpp"
#include "../../types.hpp"
#include <cmath>
#include <ostream>
#include <sstream>
#include <string>

namespace globe::io::snapshot {

// The rectangle of periods as SVG, cells drawn with straight edges. The
// frame cuts straight across the cells -- partial cells are part of the
// look -- so every cell is drawn at each period translate whose bounding
// box touches the frame and the frame clips the lot.
class FlatSvgWriter {
 public:
    explicit FlatSvgWriter(double size = 512.0);

    void write(const Snapshot& snapshot, std::ostream& stream) const;
    [[nodiscard]] std::string to_string(const Snapshot& snapshot) const;

 private:
    double _size;

    void write_cell(const Snapshot& snapshot, const Snapshot::Cell& cell, std::ostream& stream) const;
    [[nodiscard]] std::string number(double value) const;
};

inline FlatSvgWriter::FlatSvgWriter(double size) :
    _size(size) {
}

inline std::string FlatSvgWriter::to_string(const Snapshot& snapshot) const {
    std::ostringstream stream;
    write(snapshot, stream);
    return stream.str();
}

inline void FlatSvgWriter::write(const Snapshot& snapshot, std::ostream& stream) const {
    double scale = _size / snapshot.width;
    double height = snapshot.height * scale;

    stream << "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 "
        << number(_size) << " " << number(height) << "\">\n";
    stream << "  <defs><clipPath id=\"frame\"><rect x=\"0\" y=\"0\" width=\""
        << number(_size) << "\" height=\"" << number(height) << "\"/></clipPath></defs>\n";
    stream << "  <rect x=\"0\" y=\"0\" width=\"" << number(_size) << "\" height=\""
        << number(height) << "\" fill=\"#f7f5f2\"/>\n";
    stream << "  <g clip-path=\"url(#frame)\">\n";

    for (const Snapshot::Cell& cell : snapshot.cells) {
        write_cell(snapshot, cell, stream);
    }

    stream << "  </g>\n</svg>\n";
}

inline void FlatSvgWriter::write_cell(
    const Snapshot& snapshot,
    const Snapshot::Cell& cell,
    std::ostream& stream
) const {
    double scale = _size / snapshot.width;

    for (int tile_x = -1; tile_x <= 1; ++tile_x) {
        for (int tile_y = -1; tile_y <= 1; ++tile_y) {
            double offset_x = tile_x * snapshot.width;
            double offset_y = tile_y * snapshot.height;
            bool visible = false;

            for (const Vector3& vertex : cell.boundary) {
                double x = vertex.x() + offset_x;
                double y = vertex.y() + offset_y;
                visible = visible ||
                    (x > -0.25 * snapshot.width && x < 1.25 * snapshot.width &&
                     y > -0.25 * snapshot.height && y < 1.25 * snapshot.height);
            }

            if (!visible) {
                continue;
            }

            std::ostringstream path;

            for (size_t index = 0; index < cell.boundary.size(); ++index) {
                double x = (cell.boundary[index].x() + offset_x) * scale;
                double y = (snapshot.height - cell.boundary[index].y() - offset_y) * scale;
                path << (index == 0 ? "M " : "L ") << number(x) << " " << number(y) << " ";
            }

            path << "Z";

            stream << "    <path d=\"" << path.str()
                << "\" fill=\"none\" stroke=\"#2a2a2a\" stroke-width=\"1.2\"/>\n";
        }
    }
}

inline std::string FlatSvgWriter::number(double value) const {
    std::ostringstream stream;
    stream.precision(2);
    stream << std::fixed << value;
    return stream.str();
}

} // namespace globe::io::snapshot

#endif //GLOBEART_SRC_GLOBE_IO_SNAPSHOT_FLAT_SVG_WRITER_HPP_
