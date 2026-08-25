#ifndef GEOMETRY_ART_IO_QT_VIEWER_HPP_
#define GEOMETRY_ART_IO_QT_VIEWER_HPP_

#include <CGAL/Basic_viewer.h>
#include <CGAL/Graphics_scene.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Qt/qglviewer.h>
#include <QtWidgets/QWidget>
#include <QtGui/QKeyEvent>
#include <functional>
#include "../../cgal/types.hpp"
#include "../../types.hpp"
#include "../../geometry/spherical/arc.hpp"

namespace geometry_art::io::qt {

using Color = CGAL::IO::Color;

const Color BLACK(0, 0, 0);
const Color BLUE(0, 0, 255);
const Color RED(255, 0, 0);

constexpr int ARC_SEGMENTS = 20;

class Viewer final : public ::CGAL::Qt::Basic_viewer {
public:
    using KeyPressCallback = std::function<void(Viewer &, QKeyEvent *)>;

    explicit Viewer(
        QWidget *parent,
        const std::string &window_title = "Viewer",
        KeyPressCallback key_press_callback = [](Viewer &, QKeyEvent *) {}
    );

    void set_key_press_callback(const KeyPressCallback &key_press_callback) {
        _key_press_callback = key_press_callback;
    }

    void add_point(const VectorS2 &point, const Color &color = BLACK);
    void add_segment(const VectorS2 &source, const VectorS2 &target, const Color &color = BLACK);
    void add_arc(const Arc &arc, const Color &color = BLACK);
    void add_triangle(const VectorS2 &p1, const VectorS2 &p2, const VectorS2 &p3, const Color &color = BLACK);
    void add_text(const VectorS2 &point, const std::string &text);
    void clear();
    void show();
    void redraw() override;
    void center_on_unit_sphere();
    void set_edge_size(float size);

protected:
    void keyPressEvent(QKeyEvent *event) override;

private:
    KeyPressCallback _key_press_callback;
    CGAL::Graphics_scene _scene;
};

inline Viewer::Viewer(
    QWidget *parent,
    const std::string &window_title,
    KeyPressCallback key_press_callback
)
    : ::CGAL::Qt::Basic_viewer(parent, _scene, window_title.c_str()),
    _key_press_callback(key_press_callback) {
    draw_vertices(true);
    size_vertices(5.0f);
    setAxisIsDrawn(false);
    setGridIsDrawn(false);
}

inline void Viewer::add_point(const VectorS2 &point, const Color &color) {
    _scene.add_point(cgal::to_point(point), color);
}

inline void Viewer::add_segment(const VectorS2 &source, const VectorS2 &target, const Color &color) {
    _scene.add_segment(cgal::to_point(source), cgal::to_point(target), color);
}

inline void Viewer::add_arc(const Arc &arc, const Color &color) {
    for (int i = 0; i < ARC_SEGMENTS; ++i) {
        double t1 = static_cast<double>(i) / ARC_SEGMENTS;
        double t2 = static_cast<double>(i + 1) / ARC_SEGMENTS;

        add_segment(arc.interpolate(t1), arc.interpolate(t2), color);
    }
}

inline void Viewer::add_triangle(
    const VectorS2 &p1,
    const VectorS2 &p2,
    const VectorS2 &p3,
    const Color &color
) {
    VectorS2 normal = (p2 - p1).cross(p3 - p1).normalized();
    cgal::Vector3 n(normal.x(), normal.y(), normal.z());

    _scene.face_begin(color);
    _scene.add_point_in_face(cgal::to_point(p1), n);
    _scene.add_point_in_face(cgal::to_point(p2), n);
    _scene.add_point_in_face(cgal::to_point(p3), n);
    _scene.face_end();
}

inline void Viewer::add_text(const VectorS2 &point, const std::string &text) {
    _scene.add_text(cgal::to_point(point), text);
}

inline void Viewer::clear() {
    _scene.clear();
}

inline void Viewer::show() {
    ::CGAL::Qt::Basic_viewer::show();
}

inline void Viewer::redraw() {
    ::CGAL::Qt::Basic_viewer::redraw();
}

inline void Viewer::center_on_unit_sphere() {
    this->camera()->setSceneRadius(1.05);
    this->camera()->setSceneCenter(CGAL::qglviewer::Vec(0, 0, 0));
    this->camera()->setUpVector(CGAL::qglviewer::Vec(0, 0, 1));
    this->camera()->setPosition(CGAL::qglviewer::Vec(0, -3, 0));
    this->camera()->lookAt(CGAL::qglviewer::Vec(0, 0, 0));
}

inline void Viewer::set_edge_size(float size) {
    size_edges(size);
}

inline void Viewer::keyPressEvent(QKeyEvent *event) {
    _key_press_callback(*this, event);
}

} // namespace geometry_art::io::qt

#endif //GEOMETRY_ART_IO_QT_VIEWER_HPP_

