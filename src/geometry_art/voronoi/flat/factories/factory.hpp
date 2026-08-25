#ifndef GEOMETRY_ART_VORONOI_FLAT_FACTORIES_FACTORY_HPP_
#define GEOMETRY_ART_VORONOI_FLAT_FACTORIES_FACTORY_HPP_

#include "../core/torus.hpp"
#include "../optimizers/capacity_constrained_optimizer.hpp"
#include "../optimizers/lloyd_optimizer.hpp"
#include "../optimizers/newton_optimizer.hpp"
#include "../../optimizer_parameters.hpp"
#include "../../../fields/flat/constant_field.hpp"
#include "../../../fields/flat/field.hpp"
#include "../../../fields/flat/scalar_field.hpp"
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

namespace geometry_art::voronoi::flat {

using fields::flat::ConstantField;
using fields::flat::ImageField;
using fields::flat::NoiseField;
using fields::flat::PiecewisePolynomialField;
using geometry_art::math::Interval;

using SnapshotCallback = std::function<void(const io::snapshot::Snapshot&)>;

class Factory {
 public:
    Factory(
        int points_count,
        std::string density_field,
        size_t lloyd_passes,
        std::string warm_start,
        size_t newton_iterations,
        CapacityConstrainedParameters optimizer_parameters,
        std::optional<unsigned int> seed,
        double width,
        double height,
        std::string image_path,
        std::string geometry_name,
        Callback callback,
        SnapshotCallback snapshot_callback,
        std::chrono::milliseconds snapshot_interval,
        double density_tolerance = 0.0
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
    std::string _warm_start;
    size_t _newton_iterations;
    CapacityConstrainedParameters _optimizer_parameters;
    std::optional<unsigned int> _seed;
    double _width;
    double _height;
    std::string _image_path;
    std::string _geometry_name;
    Callback _callback;
    SnapshotCallback _snapshot_callback;
    std::chrono::milliseconds _snapshot_interval;
    double _density_tolerance;
    io::snapshot::Snapshot _snapshot;

    template<fields::flat::Field FieldType>
    [[nodiscard]] std::unique_ptr<Torus> build_with(const FieldType& field);

    template<fields::flat::Field FieldType>
    [[nodiscard]] std::unique_ptr<Torus> warm_started(
        std::unique_ptr<Torus> torus,
        const FieldType& field,
        const Callback& callback
    ) const;

    template<fields::flat::ScalarField ScalarFieldType>
    [[nodiscard]] PiecewisePolynomialField sample_to_tolerance(
        ScalarFieldType& scalar_field,
        const std::string& description
    ) const;

    template<fields::flat::ScalarField ScalarFieldType>
    [[nodiscard]] double relative_residual(
        const PiecewisePolynomialField& field,
        ScalarFieldType& scalar_field
    ) const;

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
    std::string warm_start,
    size_t newton_iterations,
    CapacityConstrainedParameters optimizer_parameters,
    std::optional<unsigned int> seed,
    double width,
    double height,
    std::string image_path,
    std::string geometry_name,
    Callback callback,
    SnapshotCallback snapshot_callback,
    std::chrono::milliseconds snapshot_interval,
    double density_tolerance
) :
    _points_count(points_count),
    _density_field(std::move(density_field)),
    _lloyd_passes(lloyd_passes),
    _warm_start(std::move(warm_start)),
    _newton_iterations(newton_iterations),
    _optimizer_parameters(optimizer_parameters),
    _seed(seed),
    _width(width),
    _height(height),
    _image_path(std::move(image_path)),
    _geometry_name(std::move(geometry_name)),
    _callback(std::move(callback)),
    _snapshot_callback(std::move(snapshot_callback)),
    _snapshot_interval(snapshot_interval),
    _density_tolerance(density_tolerance) {
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
    torus = warm_started(std::move(torus), field, callback);

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

    // The clock is stamped when a capture finishes rather than when it starts,
    // so a capture that costs more than the interval leaves the solver the
    // interval to run in instead of repeating back to back.
    auto last = std::make_shared<std::chrono::steady_clock::time_point>(std::chrono::steady_clock::now());

    return [this, &field, last](const Torus& torus) {
        _callback(torus);
        auto now = std::chrono::steady_clock::now();

        if (now - *last < _snapshot_interval) {
            return;
        }

        _snapshot_callback(io::snapshot::capture_flat(torus, field, _geometry_name));
        *last = std::chrono::steady_clock::now();
    };
}

inline std::unique_ptr<Torus> Factory::initial_torus() const {
    std::mt19937 engine(_seed.has_value() ? *_seed : std::random_device{}());
    std::uniform_real_distribution<double> across(0.0, _width);
    std::uniform_real_distribution<double> along(0.0, _height);

    std::vector<Vector2> sites;
    sites.reserve(static_cast<size_t>(_points_count));

    for (int k = 0; k < _points_count; ++k) {
        sites.emplace_back(across(engine), along(engine));
    }

    return std::make_unique<Torus>(_width, _height, std::move(sites));
}

template<fields::flat::Field FieldType>
std::unique_ptr<Torus> Factory::warm_started(
    std::unique_ptr<Torus> torus,
    const FieldType& field,
    const Callback& callback
) const {
    if (_warm_start == "newton" && _newton_iterations > 0) {
        NewtonParameters parameters = _optimizer_parameters.newton;
        parameters.max_iterations = _newton_iterations;
        NewtonOptimizer<FieldType> newton(std::move(torus), field, parameters, callback);
        torus = newton.optimize();

        const NewtonReport& report = newton.report();
        std::cout << "  " << std::setw(8) << std::left << "Newton" << std::right <<
            std::setw(3) << report.iterations << " steps:  gradient norm " <<
            std::scientific << std::setprecision(3) << report.gradient_norm <<
            ", energy " << report.cvt_energy << std::defaultfloat << std::endl;

        return torus;
    }

    if (_lloyd_passes == 0) {
        return torus;
    }

    LloydOptimizer<FieldType> lloyd(std::move(torus), field, _lloyd_passes, callback);
    torus = lloyd.optimize();

    std::cout << "  " << std::setw(8) << std::left << "Lloyd" << std::right <<
        std::setw(3) << _lloyd_passes << " passes: centroid deviation " <<
        std::scientific << std::setprecision(3) << lloyd.final_deviation() <<
        std::defaultfloat << std::endl;

    return torus;
}

inline PiecewisePolynomialField Factory::sample_noise_field() const {
    NoiseField noise(_width, _height, Interval(DENSITY_FLOOR, 1.0), _seed.has_value() ? static_cast<int>(*_seed) : 1546);
    return sample_to_tolerance(noise, "periodic noise");
}

inline PiecewisePolynomialField Factory::sample_image_field() const {
    ImageField image = ImageField::load(_image_path, Interval(DENSITY_FLOOR, 1.0), _width, _height);
    return sample_to_tolerance(image, "image");
}

// Accuracy as the input, as on the sphere: the grid starts at the cell
// scale and, when a tolerance is requested, doubles until the measured
// representation error meets it. The piecewise sampler has no telescoped
// residual to read off, so the error is measured by probing the scalar
// away from the fitting nodes.
template<fields::flat::ScalarField ScalarFieldType>
PiecewisePolynomialField Factory::sample_to_tolerance(
    ScalarFieldType& scalar_field,
    const std::string& description
) const {
    int columns = density_columns(_points_count, _width, _height);
    int rows = density_rows(_points_count, _width, _height);

    for (;;) {
        auto field = PiecewisePolynomialField::sample(_width, _height, columns, rows, scalar_field);
        double residual = relative_residual(field, scalar_field);

        std::cout << "Sampled " << description << " onto " << field.piece_count() <<
            " quadratic pieces, lowest sample " << field.lowest_sampled_value() <<
            ", residual RMS " << residual << " relative" << std::endl;

        bool refinable = columns < MAXIMUM_GRID_LINES || rows < MAXIMUM_GRID_LINES;

        if (_density_tolerance <= 0.0 || residual <= _density_tolerance || !refinable) {
            if (_density_tolerance > 0.0 && residual > _density_tolerance) {
                std::cout << "  WARNING: the requested density tolerance " << _density_tolerance <<
                    " was not met at the finest grid" << std::endl;
            }

            return field;
        }

        columns = std::min(2 * columns, MAXIMUM_GRID_LINES);
        rows = std::min(2 * rows, MAXIMUM_GRID_LINES);
    }
}

template<fields::flat::ScalarField ScalarFieldType>
double Factory::relative_residual(
    const PiecewisePolynomialField& field,
    ScalarFieldType& scalar_field
) const {
    size_t probes = std::min<size_t>(4 * field.piece_count(), 20000);
    double squared_error = 0.0;
    double squared_target = 0.0;

    for (size_t k = 0; k < probes; ++k) {
        double x = std::fmod(0.5 + 0.7548776662466927 * static_cast<double>(k), 1.0) * _width;
        double y = std::fmod(0.25 + 0.5698402909980532 * static_cast<double>(k), 1.0) * _height;
        double target = scalar_field.value(Vector2(x, y));
        double error = field.value(Vector2(x, y)) - target;

        squared_error += error * error;
        squared_target += target * target;
    }

    return squared_target > 0.0 ? std::sqrt(squared_error / squared_target) : 0.0;
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

} // namespace geometry_art::voronoi::flat

#endif //GEOMETRY_ART_VORONOI_FLAT_FACTORIES_FACTORY_HPP_
