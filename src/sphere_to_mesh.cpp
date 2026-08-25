#include "geometry_art/io/text/sphere_repository.hpp"
#include "geometry_art/io/mesh/exporter.hpp"
#include <CLI/CLI.hpp>
#include <filesystem>
#include <iostream>
#include <string>

using geometry_art::io::text::SphereRepository;
using geometry_art::io::mesh::Exporter;
using geometry_art::io::mesh::Format;

double parse_units(const std::string& units_str) {
    if (units_str == "mm" || units_str == "millimeters") {
        return 1.0;
    } else if (units_str == "m" || units_str == "meters") {
        return 1000.0;
    } else if (units_str == "in" || units_str == "inches") {
        return 25.4;
    }
    return 0.0;
}

int main(int argc, char* argv[]) {
    CLI::App app{"Convert sphere text file to 3D mesh"};

    std::string input_file;
    std::string format_str = "stl";
    std::string units_str = "mm";
    double diameter = 100.0;
    double arc_thickness = 1.0;
    double max_edge_length = 0.2;

    app.add_option("input", input_file)
        ->description("Input sphere text file")
        ->required();

    app.add_option("--format,-f", format_str)
        ->description("Output format: stl or ply")
        ->default_val("stl");

    app.add_option("--units,-u", units_str)
        ->description("Output units: mm, m, or in")
        ->default_val("mm");

    app.add_option("--diameter,-d", diameter)
        ->description("Target diameter in output units")
        ->default_val(100.0);

    app.add_option("--thickness,-t", arc_thickness)
        ->description("Arc bar thickness in output units")
        ->default_val(1.0);

    app.add_option("--resolution,-r", max_edge_length)
        ->description("Maximum edge length in output units")
        ->default_val(0.2);

    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError& e) {
        return app.exit(e);
    }

    Format format;
    std::string extension;

    if (format_str == "stl") {
        format = Format::STL;
        extension = ".stl";
    } else if (format_str == "ply") {
        format = Format::PLY;
        extension = ".ply";
    } else {
        std::cerr << "Unknown format: " << format_str << " (use 'stl' or 'ply')" << std::endl;
        return 1;
    }

    double mm_per_unit = parse_units(units_str);
    if (mm_per_unit == 0.0) {
        std::cerr << "Unknown units: " << units_str << " (use 'mm', 'm', or 'in')" << std::endl;
        return 1;
    }

    double scale = diameter / 2.0;

    double model_thickness = arc_thickness / scale;
    double model_resolution = max_edge_length / scale;

    std::filesystem::path input_path(input_file);

    if (!std::filesystem::exists(input_path)) {
        std::cerr << "File not found: " << input_path << std::endl;
        return 1;
    }

    auto sphere = SphereRepository::load(input_path.string());

    std::filesystem::path output_path = input_path;
    output_path.replace_extension(extension);

    std::cout <<
        "Configuration:" << std::endl <<
        "  Input: " << input_path << std::endl <<
        "  Output: " << output_path << std::endl <<
        "  Format: " << format_str << std::endl <<
        "  Units: " << units_str << std::endl <<
        "  Diameter: " << diameter << " " << units_str << std::endl <<
        "  Thickness: " << arc_thickness << " " << units_str << std::endl <<
        "  Resolution: " << max_edge_length << " " << units_str << std::endl <<
        "  Scale factor: " << scale << std::endl <<
        std::endl;

    Exporter::save(*sphere, output_path.string(), format, model_thickness, model_resolution, scale);

    std::cout << "Saved: " << output_path << std::endl;

    return 0;
}
