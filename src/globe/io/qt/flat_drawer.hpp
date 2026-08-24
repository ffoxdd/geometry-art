#ifndef GLOBEART_SRC_GLOBE_IO_QT_FLAT_DRAWER_HPP_
#define GLOBEART_SRC_GLOBE_IO_QT_FLAT_DRAWER_HPP_

#include "viewer.hpp"
#include "voronoi_sphere_drawer.hpp"
#include "../../types.hpp"
#include "../../voronoi/flat/core/torus.hpp"
#include <algorithm>
#include <cmath>
#include <optional>
#include <utility>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace globe::io::qt {

using voronoi::flat::Torus;

enum class FlatEmbedding {
    Cylinder,
    Torus
};

// Draws the flat tessellation rolled into space: a cylinder whose
// circumference is the domain's width and whose rim cuts straight across
// the cells -- partial cells are part of the look -- or the full doughnut
// when nothing is to be cut. Chart segments are subdivided so they follow
// the curved surface.
class FlatDrawer {
 public:
    FlatDrawer(
        const std::string& window_title = "Cylinder",
        RenderMode render_mode = RenderMode::Wireframe,
        FlatEmbedding embedding = FlatEmbedding::Cylinder
    );

    void show();
    void show(const Torus& torus);
    void update(const Torus& torus);

 private:
    static constexpr int SEGMENT_SAMPLES = 12;
    static constexpr int RIM_SAMPLES = 128;

    std::string _window_title;
    RenderMode _render_mode;
    FlatEmbedding _embedding;
    std::unique_ptr<Viewer> _viewer;
    bool _has_centered = false;

    void ensure_viewer();
    void draw(const Torus& torus);
    void draw_cell_boundaries(const Torus& torus, bool subdivided);
    void draw_faces(const Torus& torus);
    void draw_rims(const Torus& torus);
    void center_if_needed(const Torus& torus);

    [[nodiscard]] Vector3 embedded(const Vector2& point, const Torus& torus) const;
    [[nodiscard]] static std::vector<std::pair<double, double>> visible_spans(
        const Vector2& source,
        const Vector2& target,
        double height,
        bool clip
    );
};

inline FlatDrawer::FlatDrawer(
    const std::string& window_title,
    RenderMode render_mode,
    FlatEmbedding embedding
) :
    _window_title(window_title),
    _render_mode(render_mode),
    _embedding(embedding) {
}

inline void FlatDrawer::ensure_viewer() {
    if (!_viewer) {
        _viewer = std::make_unique<Viewer>(nullptr, _window_title);
    }
}

inline void FlatDrawer::show() {
    ensure_viewer();
    _viewer->show();
}

inline void FlatDrawer::show(const Torus& torus) {
    ensure_viewer();
    _viewer->clear();
    draw(torus);
    center_if_needed(torus);
    _viewer->show();
}

inline void FlatDrawer::update(const Torus& torus) {
    ensure_viewer();
    _viewer->clear();
    draw(torus);
    center_if_needed(torus);
    _viewer->redraw();
}

inline void FlatDrawer::draw(const Torus& torus) {
    switch (_render_mode) {
        case RenderMode::Wireframe:
            draw_cell_boundaries(torus, true);
            break;
        case RenderMode::Solid:
            draw_faces(torus);
            draw_cell_boundaries(torus, true);
            break;
        case RenderMode::Minimal:
            _viewer->set_edge_size(1.0f);
            draw_cell_boundaries(torus, false);
            break;
    }

    if (_embedding == FlatEmbedding::Cylinder) {
        draw_rims(torus);
    }
}

inline void FlatDrawer::draw_cell_boundaries(const Torus& torus, bool subdivided) {
    Color edge_color(60, 60, 60);
    bool clip = _embedding == FlatEmbedding::Cylinder;
    int samples = subdivided ? SEGMENT_SAMPLES : 1;

    // On the cylinder the seam is still periodic -- only the frame cuts --
    // so a bisector protruding past one rim reappears from the other, and
    // each period image is clipped to the band separately.
    for (size_t index = 0; index < torus.size(); ++index) {
        for (const auto& edge : torus.cell_edges(index)) {
            // Each bisector is shared; draw it from the smaller side only.
            if (edge.neighbor_index < index) {
                continue;
            }

            for (int tile : {-1, 0, 1}) {
                if (!clip && tile != 0) {
                    continue;
                }

                Vector2 offset(0.0, tile * torus.height());
                Vector2 source = edge.boundary.source() + offset;
                Vector2 target = edge.boundary.target() + offset;

                for (const auto& [low, high] : visible_spans(source, target, torus.height(), clip)) {
                    for (int sample = 0; sample < samples; ++sample) {
                        double from = low + (high - low) * sample / samples;
                        double to = low + (high - low) * (sample + 1) / samples;

                        _viewer->add_segment(
                            embedded(source + from * (target - source), torus),
                            embedded(source + to * (target - source), torus),
                            edge_color
                        );
                    }
                }
            }
        }
    }
}

inline void FlatDrawer::draw_faces(const Torus& torus) {
    Color face_color(200, 200, 200);
    bool clip = _embedding == FlatEmbedding::Cylinder;

    for (size_t index = 0; index < torus.size(); ++index) {
        auto cell = torus.cell(index);

        for (int tile : {-1, 0, 1}) {
            if (!clip && tile != 0) {
                continue;
            }

            std::vector<Vector2> shifted;
            shifted.reserve(cell.size());

            for (const Vector2& vertex : cell.vertices()) {
                shifted.push_back(vertex + Vector2(0.0, tile * torus.height()));
            }

            std::optional<geometry::planar::Polygon> visible(geometry::planar::Polygon(std::move(shifted)));

            if (clip) {
                visible = visible->clipped_by(Vector2(0.0, 1.0), Vector2(0.0, 0.0));

                if (visible.has_value()) {
                    visible = visible->clipped_by(Vector2(0.0, -1.0), Vector2(0.0, torus.height()));
                }
            }

            if (!visible.has_value() || visible->size() < 3) {
                continue;
            }

            const auto& vertices = visible->vertices();
            Vector2 centroid = Vector2::Zero();

            for (const Vector2& vertex : vertices) {
                centroid += vertex;
            }

            centroid /= static_cast<double>(vertices.size());

            for (size_t corner = 0; corner < vertices.size(); ++corner) {
                const Vector2& from = vertices[corner];
                const Vector2& to = vertices[(corner + 1) % vertices.size()];

                for (int sample = 0; sample < SEGMENT_SAMPLES; ++sample) {
                    double first = static_cast<double>(sample) / SEGMENT_SAMPLES;
                    double second = static_cast<double>(sample + 1) / SEGMENT_SAMPLES;

                    _viewer->add_triangle(
                        embedded(centroid, torus),
                        embedded(from + first * (to - from), torus),
                        embedded(from + second * (to - from), torus),
                        face_color
                    );
                }
            }
        }
    }
}

inline void FlatDrawer::draw_rims(const Torus& torus) {
    Color rim_color(140, 60, 40);

    for (double y : {0.0, torus.height()}) {
        for (int sample = 0; sample < RIM_SAMPLES; ++sample) {
            double from = torus.width() * sample / RIM_SAMPLES;
            double to = torus.width() * (sample + 1) / RIM_SAMPLES;

            _viewer->add_segment(
                embedded(Vector2(from, y), torus),
                embedded(Vector2(to, y), torus),
                rim_color
            );
        }
    }
}

inline Vector3 FlatDrawer::embedded(const Vector2& point, const Torus& torus) const {
    double around = TWO_PI * point.x() / torus.width();
    double radius = torus.width() / TWO_PI;

    if (_embedding == FlatEmbedding::Cylinder) {
        return Vector3(
            radius * std::cos(around),
            radius * std::sin(around),
            point.y() - 0.5 * torus.height()
        );
    }

    double along = TWO_PI * point.y() / torus.height();
    double tube = std::min(torus.height() / TWO_PI, 0.75 * radius);
    double ring = radius + tube * std::cos(along);

    return Vector3(ring * std::cos(around), ring * std::sin(around), tube * std::sin(along));
}

// The parameter spans of the segment that survive the frame's cut: the
// whole segment on a torus, the part with y inside [0, height] on a
// cylinder, where a protruding cell is sliced mid-cell.
inline std::vector<std::pair<double, double>> FlatDrawer::visible_spans(
    const Vector2& source,
    const Vector2& target,
    double height,
    bool clip
) {
    if (!clip) {
        return {{0.0, 1.0}};
    }

    double low = 0.0;
    double high = 1.0;
    double from_y = source.y();
    double to_y = target.y();
    double span = to_y - from_y;

    for (double bound : {0.0, height}) {
        bool from_inside = bound == 0.0 ? from_y >= bound : from_y <= bound;
        bool to_inside = bound == 0.0 ? to_y >= bound : to_y <= bound;

        if (from_inside && to_inside) {
            continue;
        }

        if (!from_inside && !to_inside) {
            return {};
        }

        double crossing = (bound - from_y) / span;

        if (from_inside) {
            high = std::min(high, crossing);
        } else {
            low = std::max(low, crossing);
        }
    }

    if (low >= high) {
        return {};
    }

    return {{low, high}};
}

inline void FlatDrawer::center_if_needed(const Torus& torus) {
    if (_has_centered) {
        return;
    }

    double radius = torus.width() / TWO_PI + torus.height();
    _viewer->camera()->setSceneRadius(radius);
    _viewer->camera()->setSceneCenter(CGAL::qglviewer::Vec(0, 0, 0));
    _viewer->camera()->setUpVector(CGAL::qglviewer::Vec(0, 0, 1));
    _viewer->camera()->setPosition(CGAL::qglviewer::Vec(0, -3.0 * radius, 0.6 * radius));
    _viewer->camera()->lookAt(CGAL::qglviewer::Vec(0, 0, 0));
    _has_centered = true;
}

} // namespace globe::io::qt

#endif //GLOBEART_SRC_GLOBE_IO_QT_FLAT_DRAWER_HPP_
