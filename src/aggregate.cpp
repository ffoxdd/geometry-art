#include "geometry_art/dla/factories/factory.hpp"
#include "geometry_art/dla/callback.hpp"
#include "geometry_art/io/snapshot/aggregate_json_writer.hpp"
#include "geometry_art/io/snapshot/aggregate_svg_writer.hpp"
#include <CLI/CLI.hpp>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <string>

using namespace geometry_art;
using dla::Aggregate;
using dla::Factory;
using dla::noop_callback;
using dla::Parameters;

struct Config {
    int particle_count;
    double particle_radius;
    double overlap;
    double spawn_margin;
    std::string return_mode;
    double kill_factor;
    std::optional<unsigned int> seed;
    std::string snapshot_path;
    double snapshot_interval = 0.0;
};

Config parse_arguments(int argc, char *argv[]);
void write_snapshot(const io::snapshot::AggregateSnapshot& snapshot, const std::string& path);

int main(int argc, char *argv[]) {
    Config config = parse_arguments(argc, argv);

    std::cout <<
        "Configuration:" << std::endl <<
        "  Particles: " << config.particle_count << std::endl <<
        "  Particle radius: " << config.particle_radius << std::endl <<
        "  Overlap: " << config.overlap << std::endl <<
        "  Spawn margin: " << config.spawn_margin << std::endl <<
        "  Return mode: " << config.return_mode << std::endl <<
        "  Seed: " << (config.seed.has_value() ? std::to_string(*config.seed) : "random") << std::endl <<
        std::endl;

    Parameters parameters;
    parameters.particle_count = static_cast<size_t>(config.particle_count);
    parameters.particle_radius = config.particle_radius;
    parameters.overlap = config.overlap;
    parameters.spawn_margin = config.spawn_margin;

    Factory factory(
        parameters,
        config.return_mode,
        config.kill_factor,
        config.seed,
        noop_callback(),
        [&](const io::snapshot::AggregateSnapshot& snapshot) { write_snapshot(snapshot, config.snapshot_path); },
        std::chrono::milliseconds(static_cast<long long>(config.snapshot_interval * 1000.0))
    );

    Aggregate<> aggregate = factory.build();

    std::cout << "Grew " << aggregate.size() << " particles, reach " << aggregate.reach() << std::endl;

    if (!config.snapshot_path.empty()) {
        write_snapshot(factory.snapshot(), config.snapshot_path);
        std::cout << "Snapshot: " << config.snapshot_path << ".json and " << config.snapshot_path << ".svg" << std::endl;
    }

    return 0;
}

// Written beside the target and renamed into place, so a reader polling the
// file never sees a partial one.
void write_snapshot(const io::snapshot::AggregateSnapshot& snapshot, const std::string& path) {
    std::filesystem::path target = std::filesystem::absolute(path);
    std::filesystem::create_directories(target.parent_path());

    {
        std::ofstream json(path + ".json.tmp");
        io::snapshot::AggregateJsonWriter().write(snapshot, json);
    }
    std::filesystem::rename(path + ".json.tmp", path + ".json");

    {
        std::ofstream drawing(path + ".svg.tmp");
        io::snapshot::AggregateSvgWriter().write(snapshot, drawing);
    }
    std::filesystem::rename(path + ".svg.tmp", path + ".svg");
}

Config parse_arguments(int argc, char *argv[]) {
    CLI::App app{"Diffusion-limited aggregation of equal spheres around a seed particle"};

    Config config;

    app.add_option("--particles,-n", config.particle_count)
        ->description("Number of particles in the finished aggregate, the seed included")
        ->default_val(2000)
        ->check(CLI::PositiveNumber);

    app.add_option("--particle-radius", config.particle_radius)
        ->description("Radius of every particle")
        ->default_val(1.0)
        ->check(CLI::PositiveNumber);

    app.add_option("--overlap", config.overlap)
        ->description("How deep a particle sticks into the one it hits, as a fraction of the radius")
        ->default_val(0.01)
        ->check(CLI::PositiveNumber);

    app.add_option("--spawn-margin", config.spawn_margin)
        ->description("Gap between the aggregate's reach and the spawn shell, in particle radii")
        ->default_val(2.0)
        ->check(CLI::PositiveNumber);

    app.add_option("--return-mode", config.return_mode)
        ->description("Far-field rule: exact harmonic return, or classic kill-and-respawn")
        ->check(CLI::IsMember({"harmonic", "kill"}))
        ->default_val("harmonic");

    app.add_option("--kill-factor", config.kill_factor)
        ->description("Kill radius as a multiple of the spawn shell radius, for --return-mode kill")
        ->default_val(3.0)
        ->check(CLI::Range(1.0, 1000.0));

    app.add_option("--seed", config.seed)
        ->description("Seed for the random walks; omit for a random seed");

    app.add_option("--snapshot", config.snapshot_path)
        ->description("Write the aggregate as JSON and SVG to this path without an extension");

    app.add_option("--snapshot-interval", config.snapshot_interval)
        ->description("Also rewrite the snapshot every this many seconds while running (0 = only at the end)")
        ->default_val(0.0)
        ->check(CLI::NonNegativeNumber);

    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
        std::exit(app.exit(e));
    }

    return config;
}
