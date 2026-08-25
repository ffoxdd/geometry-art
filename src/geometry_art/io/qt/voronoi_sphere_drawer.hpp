#ifndef GEOMETRY_ART_IO_QT_VORONOI_SPHERE_DRAWER_HPP_
#define GEOMETRY_ART_IO_QT_VORONOI_SPHERE_DRAWER_HPP_

#include "viewer.hpp"
#include "../../voronoi/spherical/core/sphere.hpp"
#include <memory>
#include <string>
#include <chrono>
#include <numeric>

namespace geometry_art::io::qt {

using voronoi::spherical::Sphere;

enum class RenderMode {
    Wireframe,
    Solid,
    Minimal
};

class SphereDrawer {
public:
    SphereDrawer(
        const std::string &window_title = "Sphere",
        RenderMode render_mode = RenderMode::Wireframe
    );

    void show();
    void show(const Sphere &sphere);
    void update(const Sphere &sphere);
    void set_render_mode(RenderMode mode);

private:
    std::string _window_title;
    RenderMode _render_mode;
    std::unique_ptr<Viewer> _viewer;
    bool _has_centered = false;

    void ensure_viewer();
    void draw(const Sphere &sphere);
    void draw_wireframe(const Sphere &sphere);
    void draw_solid(const Sphere &sphere);
    void draw_minimal(const Sphere &sphere);
    void center_if_needed();
};

inline SphereDrawer::SphereDrawer(
    const std::string &window_title,
    RenderMode render_mode
) :
    _window_title(window_title),
    _render_mode(render_mode) {
}

inline void SphereDrawer::set_render_mode(RenderMode mode) {
    _render_mode = mode;
}

inline void SphereDrawer::ensure_viewer() {
    if (!_viewer) {
        _viewer = std::make_unique<Viewer>(nullptr, _window_title);
    }
}

inline void SphereDrawer::show() {
    ensure_viewer();
    _viewer->show();
}

inline void SphereDrawer::show(const Sphere &sphere) {
    ensure_viewer();
    _viewer->clear();
    draw(sphere);
    _viewer->show();
}

inline void SphereDrawer::draw(const Sphere &sphere) {
    switch (_render_mode) {
        case RenderMode::Wireframe:
            draw_wireframe(sphere);
            break;
        case RenderMode::Solid:
            draw_solid(sphere);
            break;
        case RenderMode::Minimal:
            draw_minimal(sphere);
            break;
    }
}

inline void SphereDrawer::draw_wireframe(const Sphere &sphere) {
    Color edge_color(60, 60, 60);

    for (const auto &arc : sphere.arcs()) {
        _viewer->add_arc(arc, edge_color);
    }
}

inline void SphereDrawer::draw_solid(const Sphere &sphere) {
    Color face_color(200, 200, 200);

    for (const auto &cell : sphere.cells()) {
        std::vector<VectorS2> points(cell.points().begin(), cell.points().end());

        if (points.size() < 3) {
            continue;
        }

        VectorS2 centroid = std::accumulate(
            points.begin(), points.end(), VectorS2::Zero().eval()
        ) / static_cast<double>(points.size());

        for (size_t i = 0; i < points.size(); ++i) {
            size_t next = (i + 1) % points.size();
            _viewer->add_triangle(centroid, points[i], points[next], face_color);
        }
    }
}

inline void SphereDrawer::draw_minimal(const Sphere &sphere) {
    _viewer->set_edge_size(1.0f);
    Color edge_color(60, 60, 60);

    for (const auto &arc : sphere.arcs()) {
        _viewer->add_segment(arc.source(), arc.target(), edge_color);
    }
}

inline void SphereDrawer::center_if_needed() {
    if (!_has_centered) {
        _viewer->center_on_unit_sphere();
        _has_centered = true;
    }
}

inline void SphereDrawer::update(const Sphere &sphere) {
    ensure_viewer();
    _viewer->clear();
    draw(sphere);
    center_if_needed();
    _viewer->redraw();
}

} // namespace geometry_art::io::qt

#endif //GEOMETRY_ART_IO_QT_VORONOI_SPHERE_DRAWER_HPP_
