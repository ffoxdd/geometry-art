#ifndef GEOMETRY_ART_SKELETON_WELDER_HPP_
#define GEOMETRY_ART_SKELETON_WELDER_HPP_

#include "../cgal/types.hpp"
#include "../io/mesh/types.hpp"
#include "../types.hpp"
#include <array>
#include <cmath>
#include <cstddef>
#include <map>
#include <optional>
#include <utility>
#include <vector>

namespace geometry_art::skeleton {

using io::mesh::SurfaceMesh;
using io::mesh::VertexIndex;

// Assembles a mesh from triangles given by position, welding vertices that
// land within a tolerance of one another so that neighbouring cells, which
// each compute their shared corners on their own, close up into one
// surface. The tolerance also merges corners too close to survive the
// output's precision, whose edges would otherwise collapse there. Positions
// are welded at the model's own scale and stored at the output's.
class Welder {
 public:
    explicit Welder(double scale);

    [[nodiscard]] VertexIndex vertex(const Vector3& position);
    void triangle(VertexIndex a, VertexIndex b, VertexIndex c);

    [[nodiscard]] size_t refused() const { return _refused; }
    [[nodiscard]] SurfaceMesh take();

 private:
    static constexpr double TOLERANCE = 1e-6;

    using Key = std::array<long long, 3>;

    struct Known {
        Vector3 position;
        VertexIndex index;
    };

    double _scale;
    SurfaceMesh _mesh;
    std::map<Key, std::vector<Known>> _known;
    size_t _refused = 0;

    [[nodiscard]] std::optional<VertexIndex> nearby(const Vector3& position) const;
    [[nodiscard]] static Key key_of(const Vector3& position);
};

inline Welder::Welder(double scale) :
    _scale(scale) {
}

inline VertexIndex Welder::vertex(const Vector3& position) {
    if (std::optional<VertexIndex> found = nearby(position)) {
        return *found;
    }

    VertexIndex index = _mesh.add_vertex(cgal::to_point(Vector3(_scale * position)));
    _known[key_of(position)].push_back(Known{position, index});

    return index;
}

inline void Welder::triangle(VertexIndex a, VertexIndex b, VertexIndex c) {
    if (a == b || b == c || a == c) {
        return;
    }

    if (_mesh.add_face(a, b, c) == SurfaceMesh::null_face()) {
        ++_refused;
    }
}

inline SurfaceMesh Welder::take() {
    _known.clear();
    return std::exchange(_mesh, SurfaceMesh());
}

// A position within the tolerance of a known one has a key at most one
// step from its key along each axis.
inline std::optional<VertexIndex> Welder::nearby(const Vector3& position) const {
    Key center = key_of(position);

    for (long long x = -1; x <= 1; ++x) {
        for (long long y = -1; y <= 1; ++y) {
            for (long long z = -1; z <= 1; ++z) {
                auto found = _known.find(Key{center[0] + x, center[1] + y, center[2] + z});

                if (found == _known.end()) {
                    continue;
                }

                for (const Known& known : found->second) {
                    if ((known.position - position).lpNorm<Eigen::Infinity>() <= TOLERANCE) {
                        return known.index;
                    }
                }
            }
        }
    }

    return std::nullopt;
}

inline Welder::Key Welder::key_of(const Vector3& position) {
    return Key{
        std::llround(position.x() / TOLERANCE),
        std::llround(position.y() / TOLERANCE),
        std::llround(position.z() / TOLERANCE)
    };
}

} // namespace geometry_art::skeleton

#endif //GEOMETRY_ART_SKELETON_WELDER_HPP_
