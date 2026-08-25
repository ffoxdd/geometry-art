#ifndef GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_FACTORIES_FACTORY_HPP_
#define GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_FACTORIES_FACTORY_HPP_

#include "../core/sphere.hpp"
#include "../core/random_builder.hpp"
#include "../core/callback.hpp"
#include "../../../io/snapshot/snapshot.hpp"
#include "../optimizers/capacity_constrained_optimizer.hpp"
#include "../optimizers/lloyd_optimizer.hpp"
#include "../optimizers/newton_optimizer/newton_optimizer.hpp"
#include "../../../fields/scalar/field.hpp"
#include "../../../fields/scalar/image_field.hpp"
#include "../../../fields/scalar/noise_field.hpp"
#include "../../../math/interval.hpp"
#include "../../../fields/spherical/field.hpp"
#include "../../../fields/spherical/piecewise_polynomial_field.hpp"
#include "../../../fields/spherical/polynomial_field.hpp"
#include "../../../fields/spherical/polynomial_field_fitter.hpp"
#include "../../../fields/spherical/powell_sabin_projection.hpp"
#include "../../../generators/cartesian/random_point_generator.hpp"
#include "../../../generators/spherical/fibonacci_point_generator.hpp"
#include "../../../generators/spherical/random_point_generator.hpp"
#include "../../../geometry/cartesian/bounding_box_sampler/uniform_bounding_box_sampler.hpp"
#include "../../../math/interval_sampler/uniform_interval_sampler.hpp"
#include "../../../geometry/spherical/triangle_mesh.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <optional>
#include <string>
#include <utility>

namespace globe::voronoi::spherical {

using fields::scalar::ImageField;
using fields::scalar::NoiseField;
using globe::Interval;
using fields::spherical::PiecewisePolynomialField;
using fields::spherical::PolynomialField;
using fields::spherical::PolynomialFieldFitter;
using fields::spherical::PowellSabinProjection;
using geometry::spherical::TriangleMesh;

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
        Callback callback,
        SnapshotCallback snapshot_callback,
        std::chrono::milliseconds snapshot_interval,
        std::string image_path = "",
        double density_tolerance = 0.0
    );

    std::unique_ptr<Sphere> build();

    // The tessellation as plain data, captured against the field the
    // build actually used, which the caller no longer has a handle on.
    [[nodiscard]] const io::snapshot::Snapshot& snapshot() const { return _snapshot; }

    [[nodiscard]] static int spline_subdivisions(int points_count);

 private:
    static constexpr int NOISE_FIT_DEGREE = 8;
    static constexpr size_t NOISE_FIT_SAMPLES = 20000;
    static constexpr int NOISE_MESH_SUBDIVISIONS = 4;
    static constexpr double NOISE_DENSITY_FLOOR = 0.2;
    static constexpr int NOISE_MESH_DEGREE = 2;
    static constexpr int SPLINE_MINIMUM_SUBDIVISIONS = 1;
    static constexpr int SPLINE_MAXIMUM_SUBDIVISIONS = 6;

    using SeededBoundingBoxSampler = globe::UniformBoundingBoxSampler<globe::UniformIntervalSampler>;
    using SeededCartesianGenerator = generators::cartesian::RandomPointGenerator<SeededBoundingBoxSampler>;
    using SeededPointGenerator = generators::spherical::RandomPointGenerator<SeededCartesianGenerator>;

    int _points_count;
    std::string _density_field;
    size_t _lloyd_passes;
    std::string _warm_start;
    size_t _newton_iterations;
    CapacityConstrainedParameters _optimizer_parameters;
    std::optional<unsigned int> _seed;
    Callback _callback;
    SnapshotCallback _snapshot_callback;
    std::chrono::milliseconds _snapshot_interval;
    std::string _image_path;
    double _density_tolerance;
    io::snapshot::Snapshot _snapshot;

    template<fields::spherical::Field FieldType>
    [[nodiscard]] std::unique_ptr<Sphere> build_with(const FieldType& field);

    [[nodiscard]] PolynomialField create_polynomial_field() const;
    [[nodiscard]] std::unique_ptr<Sphere> build_initial() const;
    [[nodiscard]] static SeededPointGenerator seeded_point_generator(unsigned int seed);

    template<fields::spherical::Field FieldType>
    [[nodiscard]] Callback snapshotting_callback(const FieldType& field) const;

    template<fields::spherical::Field FieldType>
    [[nodiscard]] std::unique_ptr<Sphere> warm_start(std::unique_ptr<Sphere> sphere, const FieldType& field, const Callback& callback) const;

    template<fields::spherical::Field FieldType>
    [[nodiscard]] std::unique_ptr<Sphere> relax_with_lloyd(std::unique_ptr<Sphere> sphere, const FieldType& field, const Callback& callback) const;

    template<fields::spherical::Field FieldType>
    [[nodiscard]] std::unique_ptr<Sphere> relax_with_newton(std::unique_ptr<Sphere> sphere, const FieldType& field, const Callback& callback) const;

    template<fields::spherical::Field FieldType>
    [[nodiscard]] std::unique_ptr<Sphere> optimize(std::unique_ptr<Sphere> sphere, const FieldType& field, const Callback& callback) const;

    [[nodiscard]] static PiecewisePolynomialField sample_noise_field();
    [[nodiscard]] PiecewisePolynomialField sample_smooth_noise_field() const;
    [[nodiscard]] PiecewisePolynomialField sample_image_field() const;
    [[nodiscard]] PiecewisePolynomialField sample_quadratic_field() const;

    template<fields::scalar::Field ScalarFieldType>
    [[nodiscard]] PiecewisePolynomialField project_density(
        ScalarFieldType& scalar_field,
        const std::string& description
    ) const;
    [[nodiscard]] static PolynomialField fit_noise_field();
};

inline Factory::Factory(
    int points_count,
    std::string density_field,
    size_t lloyd_passes,
    std::string warm_start,
    size_t newton_iterations,
    CapacityConstrainedParameters optimizer_parameters,
    std::optional<unsigned int> seed,
    Callback callback,
    SnapshotCallback snapshot_callback,
    std::chrono::milliseconds snapshot_interval,
    std::string image_path,
    double density_tolerance
) :
    _points_count(points_count),
    _density_field(std::move(density_field)),
    _lloyd_passes(lloyd_passes),
    _warm_start(std::move(warm_start)),
    _newton_iterations(newton_iterations),
    _optimizer_parameters(optimizer_parameters),
    _seed(seed),
    _callback(std::move(callback)),
    _snapshot_callback(std::move(snapshot_callback)),
    _snapshot_interval(snapshot_interval),
    _image_path(std::move(image_path)),
    _density_tolerance(density_tolerance) {
}

inline std::unique_ptr<Sphere> Factory::build() {
    if (_density_field == "noise") {
        return build_with(sample_noise_field());
    }

    if (_density_field == "noise-smooth") {
        return build_with(sample_smooth_noise_field());
    }

    if (_density_field == "image") {
        return build_with(sample_image_field());
    }

    if (_density_field == "quadratic-piecewise") {
        return build_with(sample_quadratic_field());
    }

    return build_with(create_polynomial_field());
}

template<fields::spherical::Field FieldType>
std::unique_ptr<Sphere> Factory::build_with(const FieldType& field) {
    std::cout << "Generating " << _points_count << " random points..." << std::flush;
    auto sphere = build_initial();
    std::cout << " done" << std::endl;
    Callback callback = snapshotting_callback(field);
    callback(*sphere);

    sphere = warm_start(std::move(sphere), field, callback);
    sphere = optimize(std::move(sphere), field, callback);
    _snapshot = io::snapshot::capture(*sphere, field);

    return sphere;
}

// The optimizers report every step; snapshots are taken on a clock so the
// cost of capturing one stays a fixed fraction of the run regardless of how
// cheap the steps are.
template<fields::spherical::Field FieldType>
Callback Factory::snapshotting_callback(const FieldType& field) const {
    if (!_snapshot_callback || _snapshot_interval.count() <= 0) {
        return _callback;
    }

    // The clock is stamped when a capture finishes rather than when it starts,
    // so a capture that costs more than the interval leaves the solver the
    // interval to run in instead of repeating back to back.
    auto last = std::make_shared<std::chrono::steady_clock::time_point>(std::chrono::steady_clock::now());

    return [this, &field, last](const Sphere& sphere) {
        _callback(sphere);
        auto now = std::chrono::steady_clock::now();

        if (now - *last < _snapshot_interval) {
            return;
        }

        _snapshot_callback(io::snapshot::capture(sphere, field));
        *last = std::chrono::steady_clock::now();
    };
}

inline PolynomialField Factory::create_polynomial_field() const {
    if (_density_field == "constant") {
        return PolynomialField::constant(1.0);
    }

    if (_density_field == "linear") {
        return PolynomialField::linear(2.0, Vector3(0.0, 0.0, 2.0));
    }

    if (_density_field == "noise-fit") {
        return fit_noise_field();
    }

    Eigen::Matrix3d quadratic = Eigen::Matrix3d::Zero();
    quadratic(2, 2) = -0.9;
    return PolynomialField::quadratic(1.0, Vector3::Zero(), quadratic);
}

inline PiecewisePolynomialField Factory::sample_noise_field() {
    NoiseField noise_field(Interval(NOISE_DENSITY_FLOOR, 1.0));
    PiecewisePolynomialField field = PiecewisePolynomialField::sample(TriangleMesh::icosphere(NOISE_MESH_SUBDIVISIONS), NOISE_MESH_DEGREE, noise_field);

    std::cout << "Sampled noise onto piecewise degree-" << NOISE_MESH_DEGREE << " mesh with " <<
        field.mesh().triangles.size() << " triangles" << std::endl;

    return field;
}

// The C1 representation: quadratic pieces on the Powell-Sabin split, whose
// gradient is continuous and whose positivity is certified by the
// coefficients rather than sampled for.
inline PiecewisePolynomialField Factory::sample_smooth_noise_field() const {
    NoiseField noise_field(Interval(NOISE_DENSITY_FLOOR, 1.0));
    return project_density(noise_field, "noise");
}

inline PiecewisePolynomialField Factory::sample_image_field() const {
    ImageField image_field = ImageField::load(_image_path, Interval(NOISE_DENSITY_FLOOR, 1.0));

    return project_density(
        image_field,
        std::to_string(image_field.width()) + "x" + std::to_string(image_field.height()) + " image"
    );
}

// The mesh starts at the cell scale by the bandwidth rule and, when a
// tolerance is requested, refines until the projection's reported error
// meets it: accuracy is the input, and no one names a level.
template<fields::scalar::Field ScalarFieldType>
PiecewisePolynomialField Factory::project_density(
    ScalarFieldType& scalar_field,
    const std::string& description
) const {
    PowellSabinProjection projection;
    int coarsest = spline_subdivisions(_points_count);

    PowellSabinProjection::Result result = _density_tolerance > 0.0
        ? projection.project_to_tolerance(coarsest, SPLINE_MAXIMUM_SUBDIVISIONS, _density_tolerance, scalar_field)
        : projection.project(TriangleMesh::icosphere(coarsest), scalar_field);

    std::cout << "Projected " << description << " onto a C1 quadratic spline with " <<
        result.field.mesh().triangles.size() << " pieces, density at least " <<
        result.lowest_coefficient << ", residual RMS " <<
        result.root_mean_square_residual << " (" <<
        result.relative_residual << " relative)" << std::endl;

    if (_density_tolerance > 0.0 && result.relative_residual > _density_tolerance) {
        std::cout << "  WARNING: the requested density tolerance " << _density_tolerance <<
            " was not met at the finest mesh" << std::endl;
    }

    if (result.least_damping < 1.0) {
        std::cout << "  the projection dipped below zero; gradients damped to " <<
            result.least_damping << " to keep the density positive" << std::endl;
    }

    return std::move(result.field);
}

// The bandwidth rule: density structure finer than a cell cannot affect the
// tessellation except through its local average, and the projection keeps
// those averages faithful at any mesh no coarser than the cells. So the
// mesh is the coarsest icosphere whose edges fit inside a cell, and no one
// has to name a level.
inline int Factory::spline_subdivisions(int points_count) {
    double cell_radius = 2.0 / std::sqrt(static_cast<double>(points_count));
    double icosahedron_edge = std::acos(1.0 / std::sqrt(5.0));
    int level = static_cast<int>(std::ceil(std::log2(icosahedron_edge / cell_radius)));

    return std::clamp(level, SPLINE_MINIMUM_SUBDIVISIONS, SPLINE_MAXIMUM_SUBDIVISIONS);
}

// The quadratic field represented on the noise mesh: the same function, so
// any difference in behaviour is the representation's, not the field's.
inline PiecewisePolynomialField Factory::sample_quadratic_field() const {
    PolynomialField quadratic = create_polynomial_field();
    return PiecewisePolynomialField::sample(TriangleMesh::icosphere(NOISE_MESH_SUBDIVISIONS), NOISE_MESH_DEGREE, quadratic);
}

inline PolynomialField Factory::fit_noise_field() {
    NoiseField noise_field(Interval(NOISE_DENSITY_FLOOR, 1.0));
    PolynomialFieldFitter<> fitter(NOISE_FIT_DEGREE, NOISE_FIT_SAMPLES, generators::spherical::FibonacciPointGenerator());
    auto fit = fitter.fit(noise_field);

    std::cout << "Fitted noise to degree " << NOISE_FIT_DEGREE <<
        " polynomial, RMS residual " << fit.root_mean_square_residual <<
        ", lowest value " << fit.lowest_sampled_value << std::endl;

    if (fit.lowest_sampled_value <= 0.0) {
        std::cout << "  WARNING: the fit is not a density -- it reaches " << fit.lowest_sampled_value <<
            ". Cells there have no mass to balance, so the tessellation will be distorted around them." << std::endl;
    }

    return fit.field;
}

inline std::unique_ptr<Sphere> Factory::build_initial() const {
    if (!_seed.has_value()) {
        return RandomBuilder<>().build(_points_count);
    }

    return RandomBuilder<SeededPointGenerator>(seeded_point_generator(*_seed)).build(_points_count);
}

inline Factory::SeededPointGenerator Factory::seeded_point_generator(unsigned int seed) {
    return SeededPointGenerator(
        SeededCartesianGenerator(SeededBoundingBoxSampler(globe::UniformIntervalSampler(seed))),
        globe::UniformSphericalBoundingBoxSampler<>()
    );
}

template<fields::spherical::Field FieldType>
std::unique_ptr<Sphere> Factory::warm_start(std::unique_ptr<Sphere> sphere, const FieldType& field, const Callback& callback) const {
    if (_warm_start == "newton") {
        return relax_with_newton(std::move(sphere), field, callback);
    }

    return relax_with_lloyd(std::move(sphere), field, callback);
}

template<fields::spherical::Field FieldType>
std::unique_ptr<Sphere> Factory::relax_with_lloyd(std::unique_ptr<Sphere> sphere, const FieldType& field, const Callback& callback) const {
    if (_lloyd_passes == 0) {
        return sphere;
    }

    LloydOptimizer<FieldType> lloyd(std::move(sphere), field, _lloyd_passes, callback);
    sphere = lloyd.optimize();

    std::cout << "  " << std::setw(8) << std::left << "Lloyd" << std::right <<
        std::setw(3) << _lloyd_passes << " passes: centroid deviation " <<
        std::scientific << std::setprecision(3) << lloyd.final_deviation() <<
        std::defaultfloat << std::endl;

    return sphere;
}

template<fields::spherical::Field FieldType>
std::unique_ptr<Sphere> Factory::relax_with_newton(std::unique_ptr<Sphere> sphere, const FieldType& field, const Callback& callback) const {
    if (_newton_iterations == 0) {
        return sphere;
    }

    NewtonParameters parameters;
    parameters.max_iterations = _newton_iterations;
    NewtonOptimizer<FieldType> newton(std::move(sphere), field, parameters, callback);
    sphere = newton.optimize();

    const NewtonReport& report = newton.report();
    std::cout << "  " << std::setw(8) << std::left << "Newton" << std::right <<
        std::setw(3) << report.iterations << " steps:  gradient norm " <<
        std::scientific << std::setprecision(3) << report.gradient_norm <<
        ", energy " << report.cvt_energy << std::defaultfloat << std::endl;

    return sphere;
}

template<fields::spherical::Field FieldType>
std::unique_ptr<Sphere> Factory::optimize(std::unique_ptr<Sphere> sphere, const FieldType& field, const Callback& callback) const {
    CapacityConstrainedOptimizer<FieldType> optimizer(std::move(sphere), field, _optimizer_parameters, callback);
    sphere = optimizer.optimize();

    const auto& report = optimizer.report();
    std::cout << "  " << std::setw(8) << std::left << "Final" << std::right <<
        (report.converged ? " converged" : (report.stalled ? " stalled" : " stopped")) <<
        " after " << report.outer_iterations << " outer / " << report.inner_iterations << " inner iterations" <<
        ", capacity RMS " << std::scientific << std::setprecision(3) << report.relative_rms_capacity_error <<
        std::defaultfloat << std::endl;

    return sphere;
}

} // namespace globe::voronoi::spherical

#endif //GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_FACTORIES_FACTORY_HPP_
