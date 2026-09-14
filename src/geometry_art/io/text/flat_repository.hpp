#ifndef GEOMETRY_ART_IO_TEXT_FLAT_REPOSITORY_HPP_
#define GEOMETRY_ART_IO_TEXT_FLAT_REPOSITORY_HPP_

#include "../../types.hpp"
#include "../../geometry/planar/domain.hpp"
#include "../../voronoi/flat/core/diagram.hpp"
#include <fstream>
#include <iomanip>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

namespace geometry_art::io::text {

using geometry::planar::Closure;
using geometry::planar::Domain;
using voronoi::flat::Diagram;

// Sites on a flat domain as text: the first line carries the rectangle's
// extent, the second how each axis closes, every further line one
// canonical site.
class FlatRepository {
 public:
    static void save(const Diagram& diagram, const std::string& path) {
        std::ofstream file(path);

        if (!file) {
            throw std::runtime_error("Failed to open file for writing: " + path);
        }

        const Domain& domain = diagram.domain();

        file << std::setprecision(17);
        file << "extent " << domain.width << " " << domain.height << "\n";
        file << "closure " << name_of(domain.across) << " " << name_of(domain.along) << "\n";

        for (size_t index = 0; index < diagram.size(); ++index) {
            Vector2 site = diagram.site(index);
            file << site.x() << " " << site.y() << "\n";
        }
    }

    static std::unique_ptr<Diagram> load(const std::string& path) {
        std::ifstream file(path);

        if (!file) {
            throw std::runtime_error("Failed to open file for reading: " + path);
        }

        std::string header;
        double width = 0.0;
        double height = 0.0;

        if (!(file >> header >> width >> height) || header != "extent") {
            throw std::runtime_error("Missing extent header in: " + path);
        }

        std::string across;
        std::string along;

        if (!(file >> header >> across >> along) || header != "closure") {
            throw std::runtime_error("Missing closure header in: " + path);
        }

        std::vector<Vector2> sites;
        double x = 0.0;
        double y = 0.0;

        while (file >> x >> y) {
            sites.emplace_back(x, y);
        }

        Domain domain{width, height, closure_of(across, path), closure_of(along, path)};
        return std::make_unique<Diagram>(domain, std::move(sites));
    }

 private:
    static std::string name_of(Closure closure) {
        return closure == Closure::wrapped ? "wrapped" : "walled";
    }

    static Closure closure_of(const std::string& name, const std::string& path) {
        if (name == "wrapped") {
            return Closure::wrapped;
        }

        if (name == "walled") {
            return Closure::walled;
        }

        throw std::runtime_error("Unknown closure '" + name + "' in: " + path);
    }
};

} // namespace geometry_art::io::text

#endif //GEOMETRY_ART_IO_TEXT_FLAT_REPOSITORY_HPP_
