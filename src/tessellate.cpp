#include "geometry_art/voronoi/flat/factories/factory.hpp"
#include "geometry_art/voronoi/spherical/factories/factory.hpp"
#include "geometry_art/voronoi/spherical/core/callback.hpp"
#include "geometry_art/io/snapshot/flat_svg_writer.hpp"
#include "geometry_art/io/snapshot/json_writer.hpp"
#include "geometry_art/io/snapshot/svg_writer.hpp"
#include "geometry_art/io/text/sphere_repository.hpp"
#include "geometry_art/io/text/flat_repository.hpp"
#include <CLI/CLI.hpp>
#include <chrono>
#include <cmath>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <optional>
#include <sstream>
#include <string>

using namespace geometry_art;
using io::text::SphereRepository;
using voronoi::CapacityConstrainedParameters;
using voronoi::spherical::Factory;
using voronoi::spherical::Callback;
using voronoi::spherical::noop_callback;

struct Config {
    std::string geometry = "sphere";
    double scale = 50.0;
    std::string units = "output";
    std::optional<double> width;
    std::optional<double> height;
    std::optional<double> diameter;
    int points_count;
    std::string density_field;
    int lloyd_passes;
    std::string warm_start;
    int newton_iterations;
    int max_outer_iterations;
    int max_inner_iterations;
    double capacity_tolerance = 1e-7;
    std::string inner_solver;
    std::string newton_curvature;
    std::optional<unsigned int> seed;
    std::string output_dir;
    std::string snapshot_path;
    double snapshot_interval = 0.0;
    std::string image_path;
    double contrast = 4.0;
    double density_tolerance = 0.0;
};

bool flat_geometry(const std::string& geometry);
std::string flat_field_objection(const Config& config);
double model_width(const Config& config);
double model_height(const Config& config);
double in_model_units(const Config& config, double length);

Config parse_arguments(int argc, char *argv[]);
void write_snapshot(const geometry_art::io::snapshot::Snapshot& snapshot, const Config& config);

int run_flat(const Config& config, int argc, char *argv[]);

int main(int argc, char *argv[]) {
    Config config = parse_arguments(argc, argv);

    if (config.geometry != "sphere") {
        return run_flat(config, argc, argv);
    }

    std::cout <<
        "Configuration:" << std::endl <<
        "  Radius: " << config.scale << std::endl <<
        "  Points: " << config.points_count << std::endl <<
        "  Density: " << config.density_field << std::endl <<
        "  Warm start: " << config.warm_start << std::endl <<
        "  Lloyd passes: " << config.lloyd_passes << std::endl <<
        "  Newton iterations: " << config.newton_iterations << std::endl <<
        "  Max outer iterations: " << config.max_outer_iterations << std::endl <<
        "  Max inner iterations: " << config.max_inner_iterations << std::endl <<
        "  Capacity tolerance: " << config.capacity_tolerance << std::endl <<
        "  Inner solver: " << config.inner_solver << std::endl <<
        "  Contrast: " << config.contrast << std::endl <<
        "  Seed: " << (config.seed.has_value() ? std::to_string(*config.seed) : "random") << std::endl <<
        std::endl;

    Callback callback = noop_callback();

    CapacityConstrainedParameters optimizer_parameters;
    optimizer_parameters.max_outer_iterations = static_cast<size_t>(config.max_outer_iterations);
    optimizer_parameters.max_inner_iterations = static_cast<size_t>(config.max_inner_iterations);
    optimizer_parameters.relative_capacity_tolerance = config.capacity_tolerance;
    optimizer_parameters.inner_solver = config.inner_solver;
    optimizer_parameters.newton.curvature = config.newton_curvature;

    Factory factory(
        config.points_count,
        config.density_field,
        static_cast<size_t>(config.lloyd_passes),
        config.warm_start,
        static_cast<size_t>(config.newton_iterations),
        optimizer_parameters,
        config.seed,
        callback,
        [&](const geometry_art::io::snapshot::Snapshot& snapshot) { write_snapshot(snapshot, config); },
        std::chrono::milliseconds(static_cast<long long>(config.snapshot_interval * 1000.0)),
        config.image_path,
        config.contrast,
        config.density_tolerance
    );

    auto sphere = factory.build();

    if (!config.snapshot_path.empty()) {
        write_snapshot(factory.snapshot(), config);
        std::cout << "Snapshot: " << config.snapshot_path << ".json and " << config.snapshot_path << ".svg" << std::endl;
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

    return 0;
}

// Written beside the target and renamed into place, so a reader polling the
// file never sees a partial one. The snapshot is in model units and carries
// the scale that turns them back into the ones the run was sized in.
void write_snapshot(const geometry_art::io::snapshot::Snapshot& snapshot, const Config& config) {
    const std::string& path = config.snapshot_path;
    std::filesystem::path target = std::filesystem::absolute(path);
    geometry_art::io::snapshot::Snapshot scaled = snapshot;
    scaled.scale = config.scale;
    std::filesystem::create_directories(target.parent_path());

    {
        std::ofstream json(path + ".json.tmp");
        geometry_art::io::snapshot::JsonWriter().write(scaled, json);
    }
    std::filesystem::rename(path + ".json.tmp", path + ".json");

    {
        std::ofstream drawing(path + ".svg.tmp");

        if (snapshot.geometry == "sphere") {
            geometry_art::io::snapshot::SvgWriter().write(snapshot, drawing);
        } else {
            geometry_art::io::snapshot::FlatSvgWriter().write(snapshot, drawing);
        }
    }
    std::filesystem::rename(path + ".svg.tmp", path + ".svg");
}

// The flat pipeline: same options, a rectangle instead of a sphere. The
// torus wraps it both ways, the cylinder walls its rims, the plane walls
// all four sides, and walled cells end at the frame.
int run_flat(const Config& config, int argc, char *argv[]) {
    std::cout <<
        "Configuration:" << std::endl <<
        "  Geometry: " << config.geometry << " (" << model_width(config) << " x " << model_height(config) <<
        " model units, " << config.scale << " per model unit)" << std::endl <<
        "  Points: " << config.points_count << std::endl <<
        "  Density: " << config.density_field << std::endl <<
        "  Contrast: " << config.contrast << std::endl <<
        "  Warm start: " << config.warm_start << std::endl <<
        "  Lloyd passes: " << config.lloyd_passes << std::endl <<
        "  Capacity tolerance: " << config.capacity_tolerance << std::endl <<
        "  Seed: " << (config.seed.has_value() ? std::to_string(*config.seed) : "random") << std::endl <<
        std::endl;

    geometry_art::voronoi::flat::Callback callback = geometry_art::voronoi::flat::noop_callback();

    geometry_art::voronoi::CapacityConstrainedParameters optimizer_parameters;
    optimizer_parameters.max_outer_iterations = static_cast<size_t>(config.max_outer_iterations);
    optimizer_parameters.max_inner_iterations = static_cast<size_t>(config.max_inner_iterations);
    optimizer_parameters.relative_capacity_tolerance = config.capacity_tolerance;
    optimizer_parameters.inner_solver = "newton";
    optimizer_parameters.newton.curvature = config.newton_curvature;

    geometry_art::voronoi::flat::Factory factory(
        config.points_count,
        config.density_field,
        static_cast<size_t>(config.lloyd_passes),
        config.warm_start,
        static_cast<size_t>(config.newton_iterations),
        optimizer_parameters,
        config.seed,
        model_width(config),
        model_height(config),
        config.image_path,
        config.contrast,
        config.geometry,
        callback,
        [&](const geometry_art::io::snapshot::Snapshot& snapshot) { write_snapshot(snapshot, config); },
        std::chrono::milliseconds(static_cast<long long>(config.snapshot_interval * 1000.0)),
        config.density_tolerance
    );

    auto diagram = factory.build();

    if (!config.snapshot_path.empty()) {
        write_snapshot(factory.snapshot(), config);
        std::cout << "Snapshot: " << config.snapshot_path << ".json and " << config.snapshot_path << ".svg" << std::endl;
    }

    std::filesystem::create_directories(config.output_dir);

    auto time = std::time(nullptr);
    std::ostringstream filename;
    filename << config.output_dir << "/" << config.geometry << "_" <<
        config.points_count << "_" <<
        config.density_field << "_" <<
        std::put_time(std::localtime(&time), "%Y%m%d_%H%M%S") <<
        ".txt";

    geometry_art::io::text::FlatRepository::save(*diagram, filename.str());
    std::cout << "Saved: " << filename.str() << std::endl;

    return 0;
}

Config parse_arguments(int argc, char *argv[]) {
    CLI::App app{"Capacity-constrained tessellation on a sphere, torus, cylinder, or plane"};

    Config config;
    app.callback([&config]() {
        if (config.density_field == "image" && config.image_path.empty()) {
            throw CLI::ValidationError("--image", "-f image needs an image file");
        }

        std::string objection = flat_field_objection(config);

        if (!objection.empty()) {
            throw CLI::ValidationError("--density-field", objection);
        }

        if (config.diameter && config.geometry != "cylinder") {
            throw CLI::ValidationError("--diameter", "a diameter sizes the cylinder");
        }

        if (config.diameter) {
            config.width = M_PI * *config.diameter;
        }
    });

    app.add_option("--points,-p", config.points_count)
        ->description("Number of points to generate")
        ->default_val(10);

    app.add_option("--geometry,-g", config.geometry)
        ->description("Domain to tessellate: the sphere, a flat torus, a cylinder with walled rims, or a walled plane")
        ->check(CLI::IsMember({"sphere", "torus", "cylinder", "plane"}))
        ->default_val("sphere");

    app.add_option("--scale", config.scale)
        ->description("Output units per model unit: the sphere's radius")
        ->default_val(50.0)
        ->check(CLI::PositiveNumber);

    app.add_option("--units", config.units)
        ->description("What --width, --height and --diameter are given in: output or model units")
        ->check(CLI::IsMember({"output", "model"}))
        ->default_val("output");

    CLI::Option* width = app.add_option("--width", config.width)
        ->description("Width of the flat domain, the cylinder's circumference (default: 2 model units)")
        ->check(CLI::PositiveNumber);

    app.add_option("--diameter", config.diameter)
        ->description("The cylinder's diameter, in place of its circumference")
        ->excludes(width)
        ->check(CLI::PositiveNumber);

    app.add_option("--height", config.height)
        ->description("Height of the flat domain (default: 1 model unit)")
        ->check(CLI::PositiveNumber);

    app.add_option("--density-field,-f", config.density_field)
        ->description("Density field type")
        ->check(CLI::IsMember({"constant", "linear", "quadratic", "quadratic-piecewise", "noise", "noise-smooth", "noise-fit", "image"}))
        ->default_val("quadratic");

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

    app.add_option("--newton-curvature", config.newton_curvature)
        ->description("Curvature model for the newton inner solver")
        ->check(CLI::IsMember({"exact", "gauss-newton", "finite-difference"}))
        ->default_val("exact");

    app.add_option("--image", config.image_path)
        ->description("Equirectangular image whose darkness is the density, for -f image")
        ->check(CLI::ExistingFile);

    app.add_option("--contrast", config.contrast)
        ->description("Densest-to-sparsest density ratio of the linear field")
        ->default_val(4.0)
        ->check(CLI::Range(1.0, 1000.0));

    app.add_option("--density-tolerance", config.density_tolerance)
        ->description("Refine the density's spline until its relative representation error is below this (0 keeps the cell-scale mesh)")
        ->default_val(0.0)
        ->check(CLI::NonNegativeNumber);

    app.add_option("--seed", config.seed)
        ->description("Seed for the initial random points; omit for a random seed");

    app.add_option("--snapshot", config.snapshot_path)
        ->description("Write the tessellation as JSON and SVG to this path without an extension");

    app.add_option("--snapshot-interval", config.snapshot_interval)
        ->description("Also rewrite the snapshot every this many seconds while running (0 = only at the end)")
        ->default_val(0.0)
        ->check(CLI::NonNegativeNumber);

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

double model_width(const Config& config) {
    return config.width ? in_model_units(config, *config.width) : 2.0;
}

double model_height(const Config& config) {
    return config.height ? in_model_units(config, *config.height) : 1.0;
}

double in_model_units(const Config& config, double length) {
    return config.units == "model" ? length : length / config.scale;
}

bool flat_geometry(const std::string& geometry) {
    return geometry != "sphere";
}

// The flat family carries the densities it can integrate exactly, and a
// gradient only where a walled axis lets it rise.
std::string flat_field_objection(const Config& config) {
    if (!flat_geometry(config.geometry)) {
        return "";
    }

    if (config.density_field == "linear" && config.geometry == "torus") {
        return "a linear density rises from bottom to top, which the torus wraps";
    }

    if (config.density_field == "constant" || config.density_field == "linear" ||
        config.density_field == "noise" || config.density_field == "image") {
        return "";
    }

    return "the " + config.density_field + " density is not available on the " + config.geometry;
}
