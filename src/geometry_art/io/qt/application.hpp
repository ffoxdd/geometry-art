#ifndef GEOMETRY_ART_IO_QT_APPLICATION_HPP_
#define GEOMETRY_ART_IO_QT_APPLICATION_HPP_

#include <CGAL/Basic_viewer.h>
#include <QtWidgets/QApplication>
#include <memory>

namespace geometry_art::io::qt {

class Application {
public:
    Application(int argc, char *argv[]);
    int run();
    void process_events();

private:
    std::unique_ptr<QApplication> _app;
};

inline Application::Application(int argc, char *argv[]) {
    ::CGAL::Qt::init_ogl_context(4, 3);
    _app = std::make_unique<QApplication>(argc, argv);
}

inline int Application::run() {
    return _app->exec();
}

inline void Application::process_events() {
    _app->processEvents();
}

} // namespace geometry_art::io::qt

#endif //GEOMETRY_ART_IO_QT_APPLICATION_HPP_
