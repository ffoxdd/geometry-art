#ifndef GEOMETRY_ART_VORONOI_SPHERICAL_CORE_RANDOM_BUILDER_HPP_
#define GEOMETRY_ART_VORONOI_SPHERICAL_CORE_RANDOM_BUILDER_HPP_

#include "sphere.hpp"
#include "../../../generators/spherical/point_generator.hpp"
#include "../../../generators/spherical/random_point_generator.hpp"
#include <memory>

namespace geometry_art::voronoi::spherical {

template<generators::spherical::PointGenerator SpherePointGeneratorType = generators::spherical::RandomPointGenerator<>>
class RandomBuilder {
 public:
    RandomBuilder(SpherePointGeneratorType point_generator = SpherePointGeneratorType());
    [[nodiscard]] std::unique_ptr<Sphere> build(int point_count);

 private:
    SpherePointGeneratorType _point_generator;
};

template<generators::spherical::PointGenerator SpherePointGeneratorType>
RandomBuilder<SpherePointGeneratorType>::RandomBuilder(SpherePointGeneratorType point_generator) :
    _point_generator(std::move(point_generator)) {
}

template<generators::spherical::PointGenerator SpherePointGeneratorType>
std::unique_ptr<Sphere> RandomBuilder<SpherePointGeneratorType>::build(int point_count) {
    auto sphere = std::make_unique<Sphere>();

    for (const auto& point : _point_generator.generate(point_count)) {
        sphere->insert(cgal::to_point(point));
    }

    return sphere;
}

} // namespace geometry_art::voronoi::spherical

#endif //GEOMETRY_ART_VORONOI_SPHERICAL_CORE_RANDOM_BUILDER_HPP_

