#ifndef GEOMETRY_ART_DLA_FACTORIES_FACTORY_HPP_
#define GEOMETRY_ART_DLA_FACTORIES_FACTORY_HPP_

#include "../aggregate.hpp"
#include "../callback.hpp"
#include "../grower.hpp"
#include "../parameters.hpp"
#include "../particle.hpp"
#include "../walker.hpp"
#include "../direction_sampler/uniform_direction_sampler.hpp"
#include "../spawn_shell/spawn_shell.hpp"
#include "../spawn_shell/harmonic_spawn_shell.hpp"
#include "../spawn_shell/killing_spawn_shell.hpp"
#include "../../io/snapshot/aggregate_snapshot.hpp"
#include "../../math/interval_sampler/uniform_interval_sampler.hpp"
#include "../../types.hpp"
#include <chrono>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <utility>

namespace geometry_art::dla {

using SnapshotCallback = std::function<void(const io::snapshot::AggregateSnapshot&)>;

class Factory {
 public:
    Factory(
        Parameters parameters,
        std::string return_mode,
        double kill_factor,
        std::optional<unsigned int> seed,
        Callback<> callback,
        SnapshotCallback snapshot_callback,
        std::chrono::milliseconds snapshot_interval
    );

    [[nodiscard]] Aggregate<> build();

    // The finished aggregate as plain data, captured before build returned it.
    [[nodiscard]] const io::snapshot::AggregateSnapshot& snapshot() const { return _snapshot; }

 private:
    Parameters _parameters;
    std::string _return_mode;
    double _kill_factor;
    std::optional<unsigned int> _seed;
    Callback<> _callback;
    SnapshotCallback _snapshot_callback;
    std::chrono::milliseconds _snapshot_interval;
    io::snapshot::AggregateSnapshot _snapshot;

    template<SpawnShell SpawnShellType>
    [[nodiscard]] Aggregate<> build_with(SpawnShellType spawn_shell);

    [[nodiscard]] Aggregate<> initial_aggregate() const;
    [[nodiscard]] Callback<> snapshotting_callback() const;
    [[nodiscard]] UniformDirectionSampler<> direction_sampler(unsigned int stream) const;
    [[nodiscard]] UniformIntervalSampler interval_sampler(unsigned int stream) const;
};

inline Factory::Factory(
    Parameters parameters,
    std::string return_mode,
    double kill_factor,
    std::optional<unsigned int> seed,
    Callback<> callback,
    SnapshotCallback snapshot_callback,
    std::chrono::milliseconds snapshot_interval
) :
    _parameters(parameters),
    _return_mode(std::move(return_mode)),
    _kill_factor(kill_factor),
    _seed(seed),
    _callback(std::move(callback)),
    _snapshot_callback(std::move(snapshot_callback)),
    _snapshot_interval(snapshot_interval) {
}

inline Aggregate<> Factory::build() {
    if (_return_mode == "kill") {
        return build_with(KillingSpawnShell<>(_kill_factor, direction_sampler(1)));
    }

    return build_with(HarmonicSpawnShell<>(direction_sampler(1), interval_sampler(2)));
}

template<SpawnShell SpawnShellType>
Aggregate<> Factory::build_with(SpawnShellType spawn_shell) {
    Walker<UniformDirectionSampler<>, SpawnShellType> walker(
        _parameters,
        direction_sampler(0),
        std::move(spawn_shell)
    );

    Grower<RTreeParticleIndex, UniformDirectionSampler<>, SpawnShellType> grower(
        _parameters,
        initial_aggregate(),
        std::move(walker),
        snapshotting_callback()
    );

    Aggregate<> aggregate = grower.grow();
    _snapshot = io::snapshot::capture(aggregate);

    return aggregate;
}

inline Aggregate<> Factory::initial_aggregate() const {
    return Aggregate<>(_parameters.particle_radius, Particle{Vector3::Zero(), std::nullopt});
}

// The grower reports every freeze; snapshots are taken on a clock, stamped
// when a capture finishes, so an expensive capture leaves the growth the
// interval to run in instead of repeating back to back.
inline Callback<> Factory::snapshotting_callback() const {
    if (!_snapshot_callback || _snapshot_interval.count() <= 0) {
        return _callback;
    }

    auto last = std::make_shared<std::chrono::steady_clock::time_point>(std::chrono::steady_clock::now());

    return [this, last](const Aggregate<>& aggregate) {
        _callback(aggregate);
        auto now = std::chrono::steady_clock::now();

        if (now - *last < _snapshot_interval) {
            return;
        }

        _snapshot_callback(io::snapshot::capture(aggregate));
        *last = std::chrono::steady_clock::now();
    };
}

// Each consumer of randomness gets its own engine, so seeded runs stay
// reproducible however the consumers interleave their draws.
inline UniformDirectionSampler<> Factory::direction_sampler(unsigned int stream) const {
    return UniformDirectionSampler<>(interval_sampler(stream));
}

inline UniformIntervalSampler Factory::interval_sampler(unsigned int stream) const {
    if (!_seed.has_value()) {
        return UniformIntervalSampler();
    }

    return UniformIntervalSampler(*_seed + stream);
}

} // namespace geometry_art::dla

#endif //GEOMETRY_ART_DLA_FACTORIES_FACTORY_HPP_
