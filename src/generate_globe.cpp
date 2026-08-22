#include "globe/voronoi/spherical/factories/factory.hpp"
#include "globe/voronoi/spherical/core/callback.hpp"
#include "globe/io/qt/application.hpp"
#include "globe/io/qt/voronoi_sphere_drawer.hpp"
#include "globe/io/snapshot/json_writer.hpp"
#include "globe/io/snapshot/svg_writer.hpp"
#include "globe/io/text/sphere_repository.hpp"
#include <CLI/CLI.hpp>
#include <chrono>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <optional>
#include <sstream>
#include <string>

using namespace globe;
using io::qt::Application;
using io::qt::SphereDrawer;
using io::text::SphereRepository;
using voronoi::spherical::CapacityConstrainedParameters;
using voronoi::spherical::Factory;
using voronoi::spherical::Sphere;
using voronoi::spherical::Callback;
using voronoi::spherical::noop_callback;

struct Config {
    int points_count;
    std::string density_field;
    bool perform_render;
    std::string render_mode;
    int lloyd_passes;
    std::string warm_start;
    int newton_iterations;
    int max_outer_iterations;
    int max_inner_iterations;
    double capacity_tolerance = 1e-7;
    std::string inner_solver;
    std::optional<unsigned int> seed;
    std::string output_dir;
    std::string snapshot_path;
};

Config parse_arguments(int argc, char *argv[]);
void write_snapshot(const globe::io::snapshot::Snapshot& snapshot, const std::string& path);

int main(int argc, char *argv[]) {
    Config config = parse_arguments(argc, argv);

    std::cout <<
        "Configuration:" << std::endl <<
        "  Points: " << config.points_count << std::endl <<
        "  Density: " << config.density_field << std::endl <<
        "  Render: " << (config.perform_render ? "yes" : "no") << std::endl <<
        "  Render mode: " << config.render_mode << std::endl <<
        "  Warm start: " << config.warm_start << std::endl <<
        "  Lloyd passes: " << config.lloyd_passes << std::endl <<
        "  Newton iterations: " << config.newton_iterations << std::endl <<
        "  Max outer iterations: " << config.max_outer_iterations << std::endl <<
        "  Max inner iterations: " << config.max_inner_iterations << std::endl <<
        "  Capacity tolerance: " << config.capacity_tolerance << std::endl <<
        "  Inner solver: " << config.inner_solver << std::endl <<
        "  Seed: " << (config.seed.has_value() ? std::to_string(*config.seed) : "random") << std::endl <<
        std::endl;

    std::unique_ptr<Application> application;
    std::unique_ptr<SphereDrawer> drawer;
    Callback callback = noop_callback();

    if (config.perform_render) {
        application = std::make_unique<Application>(argc, argv);

        io::qt::RenderMode render_mode = io::qt::RenderMode::Wireframe;
        if (config.render_mode == "solid") {
            render_mode = io::qt::RenderMode::Solid;
        } else if (config.render_mode == "minimal") {
            render_mode = io::qt::RenderMode::Minimal;
        }

        drawer = std::make_unique<SphereDrawer>("Globe", render_mode);
        drawer->show();
        application->process_events();

        using namespace std::chrono_literals;
        auto last_render = std::chrono::steady_clock::now();

        callback = [&](const Sphere& sphere) {
            application->process_events();

            auto now = std::chrono::steady_clock::now();

            if (now - last_render >= 100ms) {
                drawer->update(sphere);
                last_render = now;
            }
        };
    }

    CapacityConstrainedParameters optimizer_parameters;
    optimizer_parameters.max_outer_iterations = static_cast<size_t>(config.max_outer_iterations);
    optimizer_parameters.max_inner_iterations = static_cast<size_t>(config.max_inner_iterations);
    optimizer_parameters.relative_capacity_tolerance = config.capacity_tolerance;
    optimizer_parameters.inner_solver = config.inner_solver;

    Factory factory(
        config.points_count,
        config.density_field,
        static_cast<size_t>(config.lloyd_passes),
        config.warm_start,
        static_cast<size_t>(config.newton_iterations),
        optimizer_parameters,
        config.seed,
        callback
    );

    auto sphere = factory.build();

    if (!config.snapshot_path.empty()) {
        write_snapshot(factory.snapshot(), config.snapshot_path);
    }

    std::filesystem::create_directories(config.output_dir);

    auto time = std::time(nullptr);
    std::ostringstream filename;
    filename << config.output_dir << "/sphere_" <<
        config.points_count << "_" <<
        config.density_field << "_" <<
        std::put_time(std::localtime(&time), "%Y%m%d_%H%M%S") <<
        ".txt";

    SphereRepository::save(*sphere, filename.str());
    std::cout << "Saved: " << filename.str() << std::endl;

    if (config.perform_render) {
        drawer->show(*sphere);
        return application->run();
    }

    return 0;
}

void write_snapshot(const globe::io::snapshot::Snapshot& snapshot, const std::string& path) {
    std::filesystem::create_directories(std::filesystem::path(path).parent_path());

    std::ofstream json(path + ".json");
    globe::io::snapshot::JsonWriter().write(snapshot, json);

    std::ofstream drawing(path + ".svg");
    globe::io::snapshot::SvgWriter().write(snapshot, drawing);

    std::cout << "Snapshot: " << path << ".json and " << path << ".svg" << std::endl;
}

Config parse_arguments(int argc, char *argv[]) {
    CLI::App app{"Globe Art Generator"};

    Config config;

    app.add_option("--points,-p", config.points_count)
        ->description("Number of points to generate")
        ->default_val(10);

    app.add_option("--density-field,-f", config.density_field)
        ->description("Density field type")
        ->check(CLI::IsMember({"constant", "linear", "quadratic", "noise", "noise-fit"}))
        ->default_val("quadratic");

    app.add_option("--render", config.perform_render)
        ->description("Enable Qt rendering")
        ->default_val(true);

    app.add_option("--render-mode", config.render_mode)
        ->description("Render mode: wireframe, solid, or minimal")
        ->check(CLI::IsMember({"wireframe", "solid", "minimal"}))
        ->default_val("wireframe");

    app.add_option("--warm-start", config.warm_start)
        ->description("Warm start method before the capacity phase")
        ->check(CLI::IsMember({"lloyd", "newton"}))
        ->default_val("lloyd");

    app.add_option("--newton-iterations", config.newton_iterations)
        ->description("Maximum trust-region Newton steps when warm starting with newton")
        ->default_val(50)
        ->check(CLI::NonNegativeNumber);

    app.add_option("--lloyd-passes", config.lloyd_passes)
        ->description("Number of density-weighted Lloyd warm-start passes")
        ->default_val(5)
        ->check(CLI::NonNegativeNumber);

    app.add_option("--max-outer-iterations", config.max_outer_iterations)
        ->description("Maximum augmented Lagrangian outer iterations")
        ->default_val(30)
        ->check(CLI::PositiveNumber);

    app.add_option("--max-inner-iterations", config.max_inner_iterations)
        ->description("Maximum L-BFGS iterations per outer iteration")
        ->default_val(200)
        ->check(CLI::PositiveNumber);

    app.add_option("--capacity-tolerance", config.capacity_tolerance)
        ->description("Relative RMS capacity error at which optimization stops")
        ->default_str("1e-7")
        ->check(CLI::PositiveNumber);

    app.add_option("--inner-solver", config.inner_solver)
        ->description("Inner minimiser for each augmented Lagrangian subproblem")
        ->check(CLI::IsMember({"lbfgs", "newton"}))
        ->default_val("lbfgs");

    app.add_option("--seed", config.seed)
        ->description("Seed for the initial random points; omit for a random seed");

    app.add_option("--snapshot", config.snapshot_path)
        ->description("Write the tessellation as JSON and SVG to this path without an extension");

    app.add_option("--output-dir,-o", config.output_dir)
        ->description("Output directory for saved spheres")
        ->default_val("./output");

    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
        std::exit(app.exit(e));
    }

    return config;
}
