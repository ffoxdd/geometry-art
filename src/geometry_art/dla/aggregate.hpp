#ifndef GEOMETRY_ART_DLA_AGGREGATE_HPP_
#define GEOMETRY_ART_DLA_AGGREGATE_HPP_

#include "particle.hpp"
#include "particle_index/particle_index.hpp"
#include "particle_index/rtree_particle_index.hpp"
#include "../types.hpp"
#include <CGAL/assertions.h>
#include <cstddef>
#include <utility>
#include <vector>

namespace geometry_art::dla {

// The frozen structure a walker diffuses toward: equal-radius spheres, each
// stuck to the one it hit, rooted at a single seed particle.
template<ParticleIndex ParticleIndexType = RTreeParticleIndex>
class Aggregate {
 public:
    Aggregate(double particle_radius, const Particle& seed);
    Aggregate(double particle_radius, const Particle& seed, ParticleIndexType index);

    [[nodiscard]] double clearance(const Vector3& position) const;
    [[nodiscard]] NearestParticle nearest_particle(const Vector3& position) const;
    void freeze(const Particle& particle);

    [[nodiscard]] const std::vector<Particle>& particles() const;
    [[nodiscard]] std::size_t size() const;
    [[nodiscard]] double particle_radius() const;
    [[nodiscard]] double contact_distance() const;
    [[nodiscard]] double reach() const;

 private:
    double _particle_radius;
    std::vector<Particle> _particles;
    ParticleIndexType _index;
    double _reach;
};

template<ParticleIndex ParticleIndexType>
Aggregate<ParticleIndexType>::Aggregate(double particle_radius, const Particle& seed) :
    Aggregate(particle_radius, seed, ParticleIndexType()) {
}

template<ParticleIndex ParticleIndexType>
Aggregate<ParticleIndexType>::Aggregate(
    double particle_radius,
    const Particle& seed,
    ParticleIndexType index
) :
    _particle_radius(particle_radius),
    _index(std::move(index)),
    _reach(0.0) {
    CGAL_precondition(particle_radius > 0.0);
    freeze(seed);
}

// The distance a walker center can travel in any direction before touching
// the aggregate.
template<ParticleIndex ParticleIndexType>
double Aggregate<ParticleIndexType>::clearance(const Vector3& position) const {
    return nearest_particle(position).distance - contact_distance();
}

template<ParticleIndex ParticleIndexType>
NearestParticle Aggregate<ParticleIndexType>::nearest_particle(const Vector3& position) const {
    return _index.nearest(position);
}

template<ParticleIndex ParticleIndexType>
void Aggregate<ParticleIndexType>::freeze(const Particle& particle) {
    _index.insert(particle.center, _particles.size());
    _particles.push_back(particle);
    _reach = std::max(_reach, particle.center.norm() + contact_distance());
}

template<ParticleIndex ParticleIndexType>
const std::vector<Particle>& Aggregate<ParticleIndexType>::particles() const {
    return _particles;
}

template<ParticleIndex ParticleIndexType>
std::size_t Aggregate<ParticleIndexType>::size() const {
    return _particles.size();
}

template<ParticleIndex ParticleIndexType>
double Aggregate<ParticleIndexType>::particle_radius() const {
    return _particle_radius;
}

// The center-to-center distance at which two particles touch.
template<ParticleIndex ParticleIndexType>
double Aggregate<ParticleIndexType>::contact_distance() const {
    return 2.0 * _particle_radius;
}

// The radius of the origin-centered sphere beyond which a walker center is
// clear of every particle.
template<ParticleIndex ParticleIndexType>
double Aggregate<ParticleIndexType>::reach() const {
    return _reach;
}

} // namespace geometry_art::dla

#endif //GEOMETRY_ART_DLA_AGGREGATE_HPP_
