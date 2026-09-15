#ifndef GEOMETRY_ART_SKELETON_WELDER_HPP_
#define GEOMETRY_ART_SKELETON_WELDER_HPP_

#include "../cgal/types.hpp"
#include "../io/mesh/types.hpp"
#include "../types.hpp"
#include <array>
#include <cmath>
#include <cstddef>
#include <map>
#include <utility>

namespace geometry_art::skeleton {

using io::mesh::SurfaceMesh;
using io::mesh::VertexIndex;

// Assembles a mesh from triangles given by position, welding vertices that
// land on the same spot so that neighbouring cells, which each compute
// their shared corners on their own, close up into one surface. Positions
// are welded at the model's own scale and stored at the output's.
class Welder {
 public:
    explicit Welder(double scale);

    [[nodiscard]] VertexIndex vertex(const Vector3& position);
    void triangle(VertexIndex a, VertexIndex b, VertexIndex c);

    [[nodiscard]] size_t refused() const { return _refused; }
    [[nodiscard]] SurfaceMesh take();

 private:
    static constexpr double QUANTUM = 1e-9;

    using Key = std::array<long long, 3>;

    double _scale;
    SurfaceMesh _mesh;
    std::map<Key, VertexIndex> _known;
    size_t _refused = 0;

    [[nodiscard]] static Key key_of(const Vector3& position);
};

inline Welder::Welder(double scale) :
    _scale(scale) {
}

inline VertexIndex Welder::vertex(const Vector3& position) {
    Key key = key_of(position);
    auto found = _known.find(key);

    if (found != _known.end()) {
        return found->second;
    }

    VertexIndex index = _mesh.add_vertex(cgal::to_point(Vector3(_scale * position)));
    _known.emplace(key, index);

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

inline Welder::Key Welder::key_of(const Vector3& position) {
    return Key{
        std::llround(position.x() / QUANTUM),
        std::llround(position.y() / QUANTUM),
        std::llround(position.z() / QUANTUM)
    };
}

} // namespace geometry_art::skeleton

#endif //GEOMETRY_ART_SKELETON_WELDER_HPP_
