#include "geometry_art/geometry/planar/domain.hpp"
#include "geometry_art/io/mesh/types.hpp"
#include "geometry_art/io/mesh/writer.hpp"
#include "geometry_art/io/snapshot/json_reader.hpp"
#include "geometry_art/io/snapshot/snapshot.hpp"
#include "geometry_art/skeleton/builder.hpp"
#include "geometry_art/skeleton/flat_outliner.hpp"
#include "geometry_art/skeleton/spherical_outliner.hpp"
#include "geometry_art/voronoi/flat/core/diagram.hpp"
#include "geometry_art/voronoi/spherical/core/sphere.hpp"
#include <CGAL/boost/graph/helpers.h>
#include <CLI/CLI.hpp>
#include <exception>
#include <filesystem>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

using namespace geometry_art;
using geometry::planar::Domain;
using io::mesh::SurfaceMesh;
using io::mesh::format_named;
using io::mesh::format_names;
using io::mesh::extension_of;
using io::snapshot::JsonReader;
using io::snapshot::Snapshot;
using skeleton::Builder;
using skeleton::FlatOutliner;
using skeleton::Parameters;
using skeleton::SphericalOutliner;
using skeleton::Window;
using voronoi::flat::Diagram;
using voronoi::spherical::Sphere;

const double DEFAULT_SCALE = 50.0;

struct Config {
    std::string snapshot_path;
    std::string output_path;
    std::string format = "stl";
    std::optional<double> scale;
    double bar_width = 1.5;
    double bar_thickness = 1.5;
    double resolution = 1.0;
    std::optional<double> window;
    std::optional<double> frame_width;
};

Config parse_arguments(int argc, char* argv[]);
std::filesystem::path output_path_of(const Config& config);
double scale_of(const Config& config, const Snapshot& snapshot);
std::optional<Window> window_of(const Config& config, double scale);
SurfaceMesh skeleton_of(const Snapshot& snapshot, const Parameters& parameters, const std::optional<Window>& window);

int main(int argc, char* argv[]) {
    Config config = parse_arguments(argc, argv);
    std::filesystem::path output = output_path_of(config);

    try {
        Snapshot snapshot = JsonReader().read_file(config.snapshot_path);
        double scale = scale_of(config, snapshot);

        Parameters parameters{
            config.bar_width / scale,
            config.bar_thickness / scale,
            config.resolution / scale,
            scale
        };

        std::optional<Window> window = window_of(config, scale);

        std::cout <<
            "Configuration:" << std::endl <<
            "  Snapshot: " << config.snapshot_path << std::endl <<
            "  Output: " << output.string() << std::endl <<
            "  Scale: " << scale << " per model unit" << std::endl <<
            "  Bars: " << config.bar_width << " wide, " << config.bar_thickness << " thick" << std::endl <<
            "  Resolution: " << config.resolution << std::endl;

        if (window) {
            std::cout << "  Window: " << *config.window << " square, framed " <<
                config.frame_width.value_or(config.bar_width) << " wide" << std::endl;
        }

        std::cout << std::endl <<
            "Skeletonizing " << snapshot.cells.size() << " cells on the " << snapshot.geometry << "..." << std::flush;

        SurfaceMesh mesh = skeleton_of(snapshot, parameters, window);
        std::cout << " done" << std::endl <<
            "  " << mesh.number_of_vertices() << " vertices, " << mesh.number_of_faces() << " faces, " <<
            (CGAL::is_closed(mesh) ? "closed" : "open") << std::endl;

        io::mesh::write(mesh, output);
        std::cout << "Saved: " << output.string() << std::endl;
    } catch (const std::exception& error) {
        std::cerr << "skeletonize: " << error.what() << std::endl;
        return 1;
    }

    return 0;
}

Config parse_arguments(int argc, char* argv[]) {
    CLI::App app{"Thicken a tessellation's edges into a printable solid"};
    Config config;

    app.add_option("snapshot", config.snapshot_path)
        ->description("The tessellation's snapshot .json")
        ->required()
        ->check(CLI::ExistingFile);

    app.add_option("--output,-o", config.output_path)
        ->description("Where to write the model; its extension picks the format (default: beside the snapshot)");

    app.add_option("--format,-f", config.format)
        ->description("Format of the model when --output is not given")
        ->check(CLI::IsMember(format_names()))
        ->default_val("stl");

    app.add_option("--scale,-s", config.scale)
        ->description("Output units per model unit: the sphere's radius, or one unit of a flat domain (default: the run's, else 50)")
        ->check(CLI::PositiveNumber);

    app.add_option("--bar-width,-w", config.bar_width)
        ->description("Width of a bar along the surface, in output units")
        ->default_val(1.5)
        ->check(CLI::PositiveNumber);

    app.add_option("--bar-thickness,-t", config.bar_thickness)
        ->description("Thickness of a bar off the surface, in output units")
        ->default_val(1.5)
        ->check(CLI::PositiveNumber);

    app.add_option("--resolution,-r", config.resolution)
        ->description("Longest facet edge along a curved surface, in output units")
        ->default_val(1.0)
        ->check(CLI::PositiveNumber);

    app.add_option("--window", config.window)
        ->description("Cut the model down to a centred square this wide, laid flat, in output units")
        ->check(CLI::PositiveNumber);

    app.add_option("--frame-width", config.frame_width)
        ->description("Width of the frame around the window, in output units (default: the bar width)")
        ->check(CLI::PositiveNumber);

    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError& error) {
        std::exit(app.exit(error));
    }

    return config;
}

std::filesystem::path output_path_of(const Config& config) {
    if (!config.output_path.empty()) {
        return config.output_path;
    }

    std::filesystem::path output = config.snapshot_path;
    output.replace_extension(extension_of(format_named(config.format)));

    return output;
}

// A run sized in output units records its scale, so a model of it comes
// out at the size it was planned at unless asked otherwise.
double scale_of(const Config& config, const Snapshot& snapshot) {
    return config.scale.value_or(snapshot.scale.value_or(DEFAULT_SCALE));
}

std::optional<Window> window_of(const Config& config, double scale) {
    if (!config.window) {
        return std::nullopt;
    }

    return Window{
        Vector2::Constant(*config.window / scale),
        config.frame_width.value_or(config.bar_width) / scale
    };
}

// The snapshot carries the sites and the domain, which is all a diagram
// needs; the cells are recomputed, since the inset needs the cuts behind
// them.
SurfaceMesh skeleton_of(const Snapshot& snapshot, const Parameters& parameters, const std::optional<Window>& window) {
    Builder builder(parameters);

    if (snapshot.geometry == "sphere" && window) {
        throw std::invalid_argument("a window needs a flat geometry");
    }

    if (snapshot.geometry == "sphere") {
        Sphere sphere;

        for (const Snapshot::Cell& cell : snapshot.cells) {
            sphere.insert(cgal::to_point(Vector3(cell.site.normalized())));
        }

        return builder.build(SphericalOutliner(sphere));
    }

    std::vector<Vector2> sites;
    sites.reserve(snapshot.cells.size());

    for (const Snapshot::Cell& cell : snapshot.cells) {
        sites.emplace_back(cell.site.x(), cell.site.y());
    }

    Diagram diagram(Domain::named(snapshot.geometry, snapshot.width, snapshot.height), std::move(sites));

    return builder.build(FlatOutliner(diagram, window));
}
