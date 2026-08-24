#ifndef GLOBEART_SRC_GLOBE_IO_TEXT_TORUS_REPOSITORY_HPP_
#define GLOBEART_SRC_GLOBE_IO_TEXT_TORUS_REPOSITORY_HPP_

#include "../../types.hpp"
#include "../../voronoi/flat/core/torus.hpp"
#include <fstream>
#include <iomanip>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

namespace globe::io::text {

using voronoi::flat::Torus;

// Sites on a flat torus as text: the first line carries the periods, every
// further line one canonical site.
class TorusRepository {
 public:
    static void save(const Torus& torus, const std::string& path) {
        std::ofstream file(path);

        if (!file) {
            throw std::runtime_error("Failed to open file for writing: " + path);
        }

        file << std::setprecision(17);
        file << "periods " << torus.width() << " " << torus.height() << "\n";

        for (size_t index = 0; index < torus.size(); ++index) {
            Vector2 site = torus.site(index);
            file << site.x() << " " << site.y() << "\n";
        }
    }

    static std::unique_ptr<Torus> load(const std::string& path) {
        std::ifstream file(path);

        if (!file) {
            throw std::runtime_error("Failed to open file for reading: " + path);
        }

        std::string header;
        double width = 0.0;
        double height = 0.0;

        if (!(file >> header >> width >> height) || header != "periods") {
            throw std::runtime_error("Missing periods header in: " + path);
        }

        auto torus = std::make_unique<Torus>(width, height);
        double x = 0.0;
        double y = 0.0;

        while (file >> x >> y) {
            torus->insert(Vector2(x, y));
        }

        return torus;
    }
};

} // namespace globe::io::text

#endif //GLOBEART_SRC_GLOBE_IO_TEXT_TORUS_REPOSITORY_HPP_
