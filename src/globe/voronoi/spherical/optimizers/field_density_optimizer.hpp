#ifndef GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_OPTIMIZERS_FIELD_DENSITY_OPTIMIZER_HPP_
#define GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_OPTIMIZERS_FIELD_DENSITY_OPTIMIZER_HPP_

#include "../../../fields/spherical/field.hpp"
#include "../../../generators/spherical/point_generator.hpp"
#include "../../../generators/spherical/random_point_generator.hpp"
#include "../../../geometry/spherical/helpers.hpp"
#include "../../../geometry/spherical/polygon/polygon.hpp"
#include "../../../types.hpp"
#include "../core/sphere.hpp"
#include "../core/callback.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <queue>
#include <utility>
#include <vector>

#include <dlib/optimization.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

namespace globe::voronoi::spherical {

using geometry::spherical::distance;

namespace detail {

struct CellMass {
    size_t index;
    double mass{};
};

struct MinMassComparator {
    bool operator()(const CellMass& a, const CellMass& b) const {
        return a.mass > b.mass;
    }
};

inline VectorS2 antipodal(const VectorS2& point) {
    return VectorS2(-point.x(), -point.y(), -point.z());
} // namespace globe::voronoi::spherical

struct TangentBasis {
    VectorS2 tangent_u;
    VectorS2 tangent_v;
};

inline bool is_pole(const VectorS2& v) {
    return std::abs(std::abs(v.z()) - 1.0) < GEOMETRIC_EPSILON;
} // namespace globe::voronoi::spherical

inline TangentBasis build_tangent_basis(const VectorS2& normal) {
    VectorS2 tangent_u;

    if (is_pole(normal)) {
        tangent_u = VectorS2(0.0, normal.z(), -normal.y());
    } else {
        tangent_u = VectorS2(normal.y(), -normal.x(), 0.0);
    }

    tangent_u = tangent_u.normalized();
    VectorS2 tangent_v = normal.cross(tangent_u);

    return {tangent_u, tangent_v};
} // namespace globe::voronoi::spherical

inline VectorS2 stereographic_plane_to_sphere(
    double u,
    double v,
    const VectorS2& south_pole,
    const VectorS2& tangent_u,
    const VectorS2& tangent_v
) {
    double r_squared = u * u + v * v;
    double scale = 4.0 / (4.0 + r_squared);

    VectorS2 result = south_pole + scale * (
        u * tangent_u + v * tangent_v - south_pole * (r_squared / 4.0)
    );

    return result.normalized();
} // namespace globe::voronoi::spherical

} // namespace globe::voronoi::spherical

template<
    fields::spherical::Field FieldType,
    generators::spherical::PointGenerator GeneratorType = generators::spherical::RandomPointGenerator<>
>
class FieldDensityOptimizer {
 public:
    static constexpr size_t DEFAULT_PASSES = 10;
    static constexpr double DISPLACEMENT_PENALTY_SCALE = 0.0;
    static constexpr double ZERO_ERROR_TOLERANCE = 1e-4;
    static constexpr double PERTURBATION_RADIANS = 0.025;
    static constexpr size_t MAX_PERTURBATION_ATTEMPTS = 100;
    static constexpr size_t MAX_SUCCESSFUL_PERTURBATIONS_BEFORE_RESTORE = 5;
    static constexpr size_t MAX_RESTORATIONS_PER_CHECKPOINT = 5;
    static constexpr double CENTROID_DEVIATION_PENALTY = 0.0;

    FieldDensityOptimizer(
        std::unique_ptr<Sphere> voronoi_sphere,
        FieldType field,
        size_t optimization_passes = DEFAULT_PASSES,
        GeneratorType point_generator = GeneratorType(),
        Callback callback = noop_callback()
    );

    std::unique_ptr<Sphere> optimize();

 private:
    std::unique_ptr<Sphere> _sphere;
    FieldType _field;
    size_t _optimization_passes;
    GeneratorType _point_generator;
    Callback _callback;

    struct PassResult {
        double start_error;
        double end_error;
        size_t vertices_processed;
    };

    enum class PassAction {
        CONTINUE,
        PERTURBING,
        RESTORED_CHECKPOINT,
        CONVERGED,
        STOPPED
    };

    struct ProgressResult {
        bool should_continue;
        PassAction action;
        size_t action_number;
    };

    struct OptimizationState {
        double best_error = std::numeric_limits<double>::max();
        std::vector<cgal::Point3> best_checkpoint;

        size_t perturbation_attempts = 0;
        size_t successful_perturbations_since_best = 0;
        size_t restorations_to_current_best = 0;

        bool just_perturbed = false;
        double error_before_perturbation = 0.0;
    };

    using CellMassHeap = std::priority_queue<
        detail::CellMass,
        std::vector<detail::CellMass>,
        detail::MinMassComparator
    >;
    using Checkpoint = std::vector<cgal::Point3>;

    double mass(const Polygon& polygon) const;
    double average_mass() const;
    double optimize_vertex_position(size_t index, double target_mass, double previous_error);
    void adjust_mass(size_t max_passes);
    PassResult run_single_pass(double target_mass);
    ProgressResult check_progress_and_maybe_perturb(
        double target_mass,
        const PassResult& pass_result,
        OptimizationState& state
    );
    void print_pass_result(
        size_t pass_number,
        size_t max_passes,
        const PassResult& result,
        const ProgressResult& progress
    );
    void print_final_results(double target_mass);
    CellMassHeap build_cell_mass_heap();
    double compute_convergence_error(double target_mass);
    double compute_objective_error(double target_mass);
    double compute_mass_error(double target_mass);
    double compute_centroid_deviation_error();
    bool perturb_most_undersized_vertex(double target_mass);
    std::optional<std::pair<size_t, double>> find_most_undersized_vertex_with_deficit(double target_mass);
    bool perturb_vertex_toward_random_point(size_t index);
    Checkpoint save_checkpoint() const;
    void restore_checkpoint(const Checkpoint& checkpoint);
};

template<fields::spherical::Field FieldType, generators::spherical::PointGenerator GeneratorType>
FieldDensityOptimizer<FieldType, GeneratorType>::FieldDensityOptimizer(
    std::unique_ptr<Sphere> voronoi_sphere,
    FieldType field,
    size_t optimization_passes,
    GeneratorType point_generator,
    Callback callback
) :
    _sphere(std::move(voronoi_sphere)),
    _field(std::move(field)),
    _optimization_passes(optimization_passes),
    _point_generator(std::move(point_generator)),
    _callback(std::move(callback)) {
}

template<fields::spherical::Field FieldType, generators::spherical::PointGenerator GeneratorType>
std::unique_ptr<Sphere> FieldDensityOptimizer<FieldType, GeneratorType>::optimize() {
    adjust_mass(_optimization_passes);
    return std::move(_sphere);
} // namespace globe::voronoi::spherical

template<fields::spherical::Field FieldType, generators::spherical::PointGenerator GeneratorType>
void FieldDensityOptimizer<FieldType, GeneratorType>::adjust_mass(size_t max_passes) {
    double target_mass = average_mass();
    std::cout << "Target mass per cell: " << target_mass << std::endl;

    OptimizationState state;
    state.best_checkpoint = save_checkpoint();

    for (size_t pass = 0; pass < max_passes; pass++) {
        PassResult result = run_single_pass(target_mass);
        ProgressResult progress = check_progress_and_maybe_perturb(target_mass, result, state);

        print_pass_result(pass + 1, max_passes, result, progress);

        if (!progress.should_continue) {
            break;
        }
    }

    print_final_results(target_mass);
} // namespace globe::voronoi::spherical

template<fields::spherical::Field FieldType, generators::spherical::PointGenerator GeneratorType>
typename FieldDensityOptimizer<FieldType, GeneratorType>::PassResult
FieldDensityOptimizer<FieldType, GeneratorType>::run_single_pass(double target_mass) {
    auto heap = build_cell_mass_heap();
    double start_error = compute_convergence_error(target_mass);
    double current_error = start_error;
    size_t vertex_count = 0;

    while (!heap.empty()) {
        size_t index = heap.top().index;
        heap.pop();

        current_error = optimize_vertex_position(index, target_mass, current_error);
        vertex_count++;
        _callback(*_sphere);
    }

    return {start_error, current_error, vertex_count};
} // namespace globe::voronoi::spherical

template<fields::spherical::Field FieldType, generators::spherical::PointGenerator GeneratorType>
typename FieldDensityOptimizer<FieldType, GeneratorType>::ProgressResult
FieldDensityOptimizer<FieldType, GeneratorType>::check_progress_and_maybe_perturb(
    double target_mass,
    const PassResult& pass_result,
    OptimizationState& state
) {
    if (state.just_perturbed) {
        state.just_perturbed = false;
        bool error_changed = std::abs(pass_result.end_error - state.error_before_perturbation) >
            ZERO_ERROR_TOLERANCE * state.error_before_perturbation;

        if (error_changed) {
            state.successful_perturbations_since_best++;
            state.perturbation_attempts = 0;
        } else {
            state.perturbation_attempts++;
        }
    }

    double pass_improvement = pass_result.start_error - pass_result.end_error;
    bool made_meaningful_progress = pass_improvement > ZERO_ERROR_TOLERANCE * pass_result.start_error;

    if (made_meaningful_progress) {
        if (pass_result.end_error < state.best_error) {
            state.best_error = pass_result.end_error;
            state.best_checkpoint = save_checkpoint();
            state.successful_perturbations_since_best = 0;
            state.restorations_to_current_best = 0;
        }
        return {true, PassAction::CONTINUE, 0};
    }

    double rms_error = std::sqrt(pass_result.end_error / _sphere->size());
    if (rms_error <= ZERO_ERROR_TOLERANCE * target_mass) {
        return {false, PassAction::CONVERGED, 0};
    }

    if (state.successful_perturbations_since_best >= MAX_SUCCESSFUL_PERTURBATIONS_BEFORE_RESTORE &&
        state.restorations_to_current_best < MAX_RESTORATIONS_PER_CHECKPOINT) {
        restore_checkpoint(state.best_checkpoint);
        state.restorations_to_current_best++;
        state.successful_perturbations_since_best = 0;
        state.perturbation_attempts = 0;
        return {true, PassAction::RESTORED_CHECKPOINT, state.restorations_to_current_best};
    }

    if (state.perturbation_attempts >= MAX_PERTURBATION_ATTEMPTS) {
        return {false, PassAction::STOPPED, state.perturbation_attempts};
    }

    state.error_before_perturbation = pass_result.end_error;
    size_t attempt_number = state.perturbation_attempts + 1;

    if (perturb_most_undersized_vertex(target_mass)) {
        state.just_perturbed = true;
        return {true, PassAction::PERTURBING, attempt_number};
    } else {
        state.perturbation_attempts++;
        return {true, PassAction::CONTINUE, 0};
    }
} // namespace globe::voronoi::spherical

template<fields::spherical::Field FieldType, generators::spherical::PointGenerator GeneratorType>
void FieldDensityOptimizer<FieldType, GeneratorType>::print_pass_result(
    size_t pass_number,
    size_t max_passes,
    const PassResult& result,
    const ProgressResult& progress
) {
    double start_rms = std::sqrt(result.start_error / _sphere->size());
    double end_rms = std::sqrt(result.end_error / _sphere->size());

    std::cout << std::fixed << std::setprecision(6) <<
        "  Pass " << std::setw(3) << pass_number << "/" << max_passes <<
        ": error " << std::setw(10) << start_rms <<
        " -> " << std::setw(10) << end_rms;

    switch (progress.action) {
        case PassAction::PERTURBING:
            std::cout << " [perturb #" << progress.action_number << "]";
            break;
        case PassAction::RESTORED_CHECKPOINT:
            std::cout << " [restore #" << progress.action_number << "]";
            break;
        case PassAction::CONVERGED:
            std::cout << " [converged]";
            break;
        case PassAction::STOPPED:
            std::cout << " [stopped]";
            break;
        case PassAction::CONTINUE:
            break;
    }

    std::cout << std::defaultfloat << std::endl;
} // namespace globe::voronoi::spherical

template<fields::spherical::Field FieldType, generators::spherical::PointGenerator GeneratorType>
void FieldDensityOptimizer<FieldType, GeneratorType>::print_final_results(double target_mass) {
    std::cout << std::endl;
    std::cout << "Final cell masses after optimization:" << std::endl;

    double total_error = 0.0;
    double max_error = 0.0;

    size_t i = 0;
    for (const auto& cell : _sphere->cells()) {
        double cell_mass = mass(cell);
        double error = std::abs(cell_mass - target_mass);

        total_error += error;
        max_error = std::max(max_error, error);

        std::cout << "  Cell " << std::setw(4) << i << " mass: " << cell_mass << ", error: " << error << std::endl;
        i++;
    }

    std::cout << "Average error: " << (total_error / _sphere->size()) << std::endl;
    std::cout << "Max error: " << max_error << std::endl;
} // namespace globe::voronoi::spherical

template<fields::spherical::Field FieldType, generators::spherical::PointGenerator GeneratorType>
double FieldDensityOptimizer<FieldType, GeneratorType>::mass(const Polygon& polygon) const {
    return _field.integrals(polygon).mass;
} // namespace globe::voronoi::spherical

template<fields::spherical::Field FieldType, generators::spherical::PointGenerator GeneratorType>
double FieldDensityOptimizer<FieldType, GeneratorType>::average_mass() const {
    return _field.total_mass() / _sphere->size();
} // namespace globe::voronoi::spherical

template<fields::spherical::Field FieldType, generators::spherical::PointGenerator GeneratorType>
double FieldDensityOptimizer<FieldType, GeneratorType>::optimize_vertex_position(
    size_t index,
    double target_mass,
    double previous_error
) {
    cgal::Point3 original_position = _sphere->site(index);
    VectorS2 original_vector = to_vector_s2(original_position);
    VectorS2 north = detail::antipodal(original_vector);
    detail::TangentBasis basis = detail::build_tangent_basis(north);
    VectorS2 tangent_u = basis.tangent_u;
    VectorS2 tangent_v = basis.tangent_v;

    using column_vector = dlib::matrix<double, 2, 1>;

    double initial_error = previous_error;

    auto run_bobyqa = [&](const column_vector& initial_point, int max_function_calls) {
        struct RunResult {
            double mass_error;
            double cost;
            cgal::Point3 best_position;
        };

        cgal::Point3 run_best_position = original_position;
        double run_best_mass_error = initial_error;
        double run_best_cost = initial_error;
        column_vector starting_point = initial_point;
        _sphere->update_site(index, original_position);

        auto objective = [&](const column_vector& params) -> double {
            VectorS2 candidate_vec = detail::stereographic_plane_to_sphere(
                params(0),
                params(1),
                original_vector,
                tangent_u,
                tangent_v
            );
            cgal::Point3 candidate = cgal::to_point(candidate_vec);

            _sphere->update_site(index, candidate);

            double optimization_error = compute_objective_error(target_mass);
            double displacement_angle = distance(original_vector, candidate_vec);
            double penalty = DISPLACEMENT_PENALTY_SCALE * target_mass * target_mass * displacement_angle * displacement_angle;
            double cost = optimization_error + penalty;

            if (cost < run_best_cost) {
                run_best_cost = cost;
                run_best_mass_error = compute_convergence_error(target_mass);
                run_best_position = candidate;
            }

            return cost;
        };

        try {
            dlib::find_min_bobyqa(
                objective,
                starting_point,
                5,
                dlib::uniform_matrix<double>(2, 1, -10.0),
                dlib::uniform_matrix<double>(2, 1, 10.0),
                0.5,
                1e-5,
                max_function_calls
            );
        } catch (...) {
        }

        _sphere->update_site(index, run_best_position);
        return RunResult{run_best_mass_error, run_best_cost, run_best_position};
    };

    column_vector origin;
    origin = 0.0, 0.0;

    auto initial_result = run_bobyqa(origin, 15);
    _sphere->update_site(index, initial_result.best_position);

    return initial_result.mass_error;
} // namespace globe::voronoi::spherical

template<fields::spherical::Field FieldType, generators::spherical::PointGenerator GeneratorType>
typename FieldDensityOptimizer<FieldType, GeneratorType>::CellMassHeap
FieldDensityOptimizer<FieldType, GeneratorType>::build_cell_mass_heap() {
    CellMassHeap heap;

    size_t i = 0;
    for (const auto& cell : _sphere->cells()) {
        double cell_mass = mass(cell);
        heap.push({i, cell_mass});
        i++;
    }

    return heap;
} // namespace globe::voronoi::spherical

template<fields::spherical::Field FieldType, generators::spherical::PointGenerator GeneratorType>
double FieldDensityOptimizer<FieldType, GeneratorType>::compute_convergence_error(double target_mass) {
    return compute_mass_error(target_mass);
} // namespace globe::voronoi::spherical

template<fields::spherical::Field FieldType, generators::spherical::PointGenerator GeneratorType>
double FieldDensityOptimizer<FieldType, GeneratorType>::compute_objective_error(double target_mass) {
    return compute_mass_error(target_mass) +
        CENTROID_DEVIATION_PENALTY * target_mass * target_mass * compute_centroid_deviation_error();
} // namespace globe::voronoi::spherical

template<fields::spherical::Field FieldType, generators::spherical::PointGenerator GeneratorType>
double FieldDensityOptimizer<FieldType, GeneratorType>::compute_mass_error(double target_mass) {
    std::vector<Polygon> cells;
    cells.reserve(_sphere->size());

    for (const auto& cell : _sphere->cells()) {
        cells.push_back(cell);
    }

    std::atomic<double> total_error{0.0};

    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, cells.size()),
        [&](const tbb::blocked_range<size_t>& range) {
            double local_error = 0.0;

            for (size_t i = range.begin(); i != range.end(); ++i) {
                double cell_mass = mass(cells[i]);
                double mass_error = cell_mass - target_mass;
                local_error += mass_error * mass_error;
            }

            double current = total_error.load();
            while (!total_error.compare_exchange_weak(current, current + local_error));
        }
    );

    return total_error.load();
} // namespace globe::voronoi::spherical

template<fields::spherical::Field FieldType, generators::spherical::PointGenerator GeneratorType>
double FieldDensityOptimizer<FieldType, GeneratorType>::compute_centroid_deviation_error() {
    double total_error = 0.0;

    size_t i = 0;
    for (const auto& cell : _sphere->cells()) {
        VectorS2 site = to_vector_s2(_sphere->site(i));
        VectorS2 centroid = cell.centroid();
        double deviation = distance(site, centroid);
        total_error += deviation * deviation;
        i++;
    }

    return total_error;
} // namespace globe::voronoi::spherical

template<fields::spherical::Field FieldType, generators::spherical::PointGenerator GeneratorType>
bool FieldDensityOptimizer<FieldType, GeneratorType>::perturb_most_undersized_vertex(
    double target_mass
) {
    auto candidate = find_most_undersized_vertex_with_deficit(target_mass);

    if (!candidate.has_value()) {
        return false;
    }

    auto [vertex_index, deficit] = candidate.value();

    if (!perturb_vertex_toward_random_point(vertex_index)) {
        return false;
    }

    return true;
} // namespace globe::voronoi::spherical

template<fields::spherical::Field FieldType, generators::spherical::PointGenerator GeneratorType>
std::optional<std::pair<size_t, double>>
FieldDensityOptimizer<FieldType, GeneratorType>::find_most_undersized_vertex_with_deficit(
    double target_mass
) {
    size_t best_index = _sphere->size();
    double largest_deficit = 0.0;
    size_t i = 0;

    for (const auto& cell : _sphere->cells()) {
        double cell_mass = mass(cell);
        double deficit = target_mass - cell_mass;

        if (deficit > largest_deficit) {
            largest_deficit = deficit;
            best_index = i;
        }

        i++;
    }

    if (largest_deficit <= 0.0) {
        return std::nullopt;
    }

    return std::make_pair(best_index, largest_deficit);
} // namespace globe::voronoi::spherical

template<fields::spherical::Field FieldType, generators::spherical::PointGenerator GeneratorType>
bool FieldDensityOptimizer<FieldType, GeneratorType>::perturb_vertex_toward_random_point(
    size_t index
) {
    if (index >= _sphere->size()) {
        return false;
    }

    VectorS2 current_vector = to_vector_s2(_sphere->site(index));
    VectorS2 random_vector = _point_generator.generate(1)[0];

    double dot = current_vector.dot(random_vector);
    VectorS2 tangent_component = random_vector - dot * current_vector;
    double tangent_length = tangent_component.norm();

    if (tangent_length < 1e-10) {
        return false;
    }

    VectorS2 tangent_direction = tangent_component / tangent_length;
    VectorS2 new_vector =
        std::cos(PERTURBATION_RADIANS) * current_vector +
        std::sin(PERTURBATION_RADIANS) * tangent_direction;

    _sphere->update_site(index, cgal::to_point(new_vector));

    return true;
} // namespace globe::voronoi::spherical

template<fields::spherical::Field FieldType, generators::spherical::PointGenerator GeneratorType>
typename FieldDensityOptimizer<FieldType, GeneratorType>::Checkpoint
FieldDensityOptimizer<FieldType, GeneratorType>::save_checkpoint() const {
    Checkpoint checkpoint;
    checkpoint.reserve(_sphere->size());

    for (size_t i = 0; i < _sphere->size(); i++) {
        checkpoint.push_back(_sphere->site(i));
    }

    return checkpoint;
} // namespace globe::voronoi::spherical

template<fields::spherical::Field FieldType, generators::spherical::PointGenerator GeneratorType>
void FieldDensityOptimizer<FieldType, GeneratorType>::restore_checkpoint(
    const Checkpoint& checkpoint
) {
    auto new_voronoi = std::make_unique<Sphere>();
    for (const auto& site : checkpoint) {
        new_voronoi->insert(site);
    }
    _sphere = std::move(new_voronoi);
} // namespace globe::voronoi::spherical

} // namespace globe::voronoi::spherical

#endif //GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_OPTIMIZERS_FIELD_DENSITY_OPTIMIZER_HPP_
