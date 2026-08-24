#ifndef GLOBEART_SRC_GLOBE_VORONOI_FLAT_FACTORIES_FACTORY_HPP_
#define GLOBEART_SRC_GLOBE_VORONOI_FLAT_FACTORIES_FACTORY_HPP_

#include "../core/torus.hpp"
#include "../optimizers/capacity_constrained_optimizer.hpp"
#include "../optimizers/lloyd_optimizer.hpp"
#include "../../optimizer_parameters.hpp"
#include "../../../fields/flat/constant_field.hpp"
#include "../../../fields/flat/field.hpp"
#include "../../../fields/flat/image_field.hpp"
#include "../../../fields/flat/noise_field.hpp"
#include "../../../fields/flat/piecewise_polynomial_field.hpp"
#include "../../../io/snapshot/flat_capture.hpp"
#include "../../../io/snapshot/snapshot.hpp"
#include "../../../math/interval.hpp"
#include "../../../types.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <optional>
#include <random>
#include <string>
#include <utility>

namespace globe::voronoi::flat {

using fields::flat::ConstantField;
using fields::flat::ImageField;
using fields::flat::NoiseField;
using fields::flat::PiecewisePolynomialField;
using globe::math::Interval;

using SnapshotCallback = std::function<void(const io::snapshot::Snapshot&)>;

class Factory {
 public:
    Factory(
        int points_count,
        std::string density_field,
        size_t lloyd_passes,
        CapacityConstrainedParameters optimizer_parameters,
        std::optional<unsigned int> seed,
        double width,
        double height,
        std::string image_path,
        std::string geometry_name,
        Callback callback,
        SnapshotCallback snapshot_callback,
        std::chrono::milliseconds snapshot_interval
    );

    std::unique_ptr<Torus> build();

    [[nodiscard]] const io::snapshot::Snapshot& snapshot() const { return _snapshot; }

    // The flat bandwidth rule: the density grid's spacing is the cell
    // scale, the square root of one cell's share of the rectangle.
    [[nodiscard]] static int density_columns(int points_count, double width, double height);
    [[nodiscard]] static int density_rows(int points_count, double width, double height);

 private:
    static constexpr double DENSITY_FLOOR = 0.2;
    static constexpr int MINIMUM_GRID_LINES = 4;
    static constexpr int MAXIMUM_GRID_LINES = 256;

    int _points_count;
    std::string _density_field;
    size_t _lloyd_passes;
    CapacityConstrainedParameters _optimizer_parameters;
    std::optional<unsigned int> _seed;
    double _width;
    double _height;
    std::string _image_path;
    std::string _geometry_name;
    Callback _callback;
    SnapshotCallback _snapshot_callback;
    std::chrono::milliseconds _snapshot_interval;
    io::snapshot::Snapshot _snapshot;

    template<fields::flat::Field FieldType>
    [[nodiscard]] std::unique_ptr<Torus> build_with(const FieldType& field);

    [[nodiscard]] std::unique_ptr<Torus> initial_torus() const;
    [[nodiscard]] PiecewisePolynomialField sample_noise_field() const;
    [[nodiscard]] PiecewisePolynomialField sample_image_field() const;

    template<fields::flat::Field FieldType>
    [[nodiscard]] Callback snapshotting_callback(const FieldType& field) const;
};

inline Factory::Factory(
    int points_count,
    std::string density_field,
    size_t lloyd_passes,
    CapacityConstrainedParameters optimizer_parameters,
    std::optional<unsigned int> seed,
    double width,
    double height,
    std::string image_path,
    std::string geometry_name,
    Callback callback,
    SnapshotCallback snapshot_callback,
    std::chrono::milliseconds snapshot_interval
) :
    _points_count(points_count),
    _density_field(std::move(density_field)),
    _lloyd_passes(lloyd_passes),
    _optimizer_parameters(optimizer_parameters),
    _seed(seed),
    _width(width),
    _height(height),
    _image_path(std::move(image_path)),
    _geometry_name(std::move(geometry_name)),
    _callback(std::move(callback)),
    _snapshot_callback(std::move(snapshot_callback)),
    _snapshot_interval(snapshot_interval) {
}

inline std::unique_ptr<Torus> Factory::build() {
    if (_density_field == "noise") {
        return build_with(sample_noise_field());
    }

    if (_density_field == "image") {
        return build_with(sample_image_field());
    }

    return build_with(ConstantField(1.0, _width, _height));
}

template<fields::flat::Field FieldType>
std::unique_ptr<Torus> Factory::build_with(const FieldType& field) {
    std::cout << "Scattering " << _points_count << " random points on the " <<
        _geometry_name << "..." << std::flush;
    auto torus = initial_torus();
    std::cout << " done" << std::endl;

    Callback callback = snapshotting_callback(field);
    callback(*torus);

    LloydOptimizer<FieldType> lloyd(std::move(torus), field, _lloyd_passes, callback);
    torus = lloyd.optimize();

    if (_lloyd_passes > 0) {
        std::cout << "  " << std::setw(8) << std::left << "Lloyd" << std::right <<
            std::setw(3) << _lloyd_passes << " passes: centroid deviation " <<
            std::scientific << std::setprecision(3) << lloyd.final_deviation() <<
            std::defaultfloat << std::endl;
    }

    CapacityConstrainedOptimizer<FieldType> optimizer(std::move(torus), field, _optimizer_parameters, callback);
    torus = optimizer.optimize();

    const auto& report = optimizer.report();
    std::cout << "  " << std::setw(8) << std::left << "Final" << std::right <<
        (report.converged ? " converged" : (report.stalled ? " stalled" : " stopped")) <<
        " after " << report.outer_iterations << " outer / " << report.inner_iterations << " inner iterations" <<
        ", capacity RMS " << std::scientific << std::setprecision(3) << report.relative_rms_capacity_error <<
        std::defaultfloat << std::endl;

    _snapshot = io::snapshot::capture_flat(*torus, field, _geometry_name);
    return torus;
}

template<fields::flat::Field FieldType>
Callback Factory::snapshotting_callback(const FieldType& field) const {
    if (!_snapshot_callback || _snapshot_interval.count() <= 0) {
        return _callback;
    }

    auto last = std::make_shared<std::chrono::steady_clock::time_point>(std::chrono::steady_clock::now());

    return [this, &field, last](const Torus& torus) {
        _callback(torus);
        auto now = std::chrono::steady_clock::now();

        if (now - *last < _snapshot_interval) {
            return;
        }

        *last = now;
        _snapshot_callback(io::snapshot::capture_flat(torus, field, _geometry_name));
    };
}

inline std::unique_ptr<Torus> Factory::initial_torus() const {
    std::mt19937 engine(_seed.has_value() ? *_seed : std::random_device{}());
    std::uniform_real_distribution<double> across(0.0, _width);
    std::uniform_real_distribution<double> along(0.0, _height);

    auto torus = std::make_unique<Torus>(_width, _height);

    for (int k = 0; k < _points_count; ++k) {
        torus->insert(Vector2(across(engine), along(engine)));
    }

    return torus;
}

inline PiecewisePolynomialField Factory::sample_noise_field() const {
    NoiseField noise(_width, _height, Interval(DENSITY_FLOOR, 1.0), _seed.has_value() ? static_cast<int>(*_seed) : 1546);
    auto field = PiecewisePolynomialField::sample(
        _width,
        _height,
        density_columns(_points_count, _width, _height),
        density_rows(_points_count, _width, _height),
        noise
    );

    std::cout << "Sampled periodic noise onto " << field.piece_count() <<
        " quadratic pieces, lowest sample " << field.lowest_sampled_value() << std::endl;

    return field;
}

inline PiecewisePolynomialField Factory::sample_image_field() const {
    ImageField image = ImageField::load(_image_path, Interval(DENSITY_FLOOR, 1.0), _width, _height);
    auto field = PiecewisePolynomialField::sample(
        _width,
        _height,
        density_columns(_points_count, _width, _height),
        density_rows(_points_count, _width, _height),
        image
    );

    std::cout << "Sampled image onto " << field.piece_count() <<
        " quadratic pieces, lowest sample " << field.lowest_sampled_value() << std::endl;

    return field;
}

inline int Factory::density_columns(int points_count, double width, double height) {
    double cell_size = std::sqrt(width * height / static_cast<double>(points_count));
    int lines = static_cast<int>(std::ceil(width / cell_size));
    return std::clamp(lines, MINIMUM_GRID_LINES, MAXIMUM_GRID_LINES);
}

inline int Factory::density_rows(int points_count, double width, double height) {
    double cell_size = std::sqrt(width * height / static_cast<double>(points_count));
    int lines = static_cast<int>(std::ceil(height / cell_size));
    return std::clamp(lines, MINIMUM_GRID_LINES, MAXIMUM_GRID_LINES);
}

} // namespace globe::voronoi::flat

#endif //GLOBEART_SRC_GLOBE_VORONOI_FLAT_FACTORIES_FACTORY_HPP_
