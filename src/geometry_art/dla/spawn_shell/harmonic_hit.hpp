#ifndef GEOMETRY_ART_DLA_SPAWN_SHELL_HARMONIC_HIT_HPP_
#define GEOMETRY_ART_DLA_SPAWN_SHELL_HARMONIC_HIT_HPP_

namespace geometry_art::dla {

// A Brownian walker at the given distance outside a sphere of the given
// radius either returns to the sphere or wanders off to infinity. Both halves
// of that fate are classical and closed-form: it returns with probability
// radius / distance, and a returning walker lands with the density of the
// exterior Poisson kernel, whose polar CDF inverts exactly. Together they
// resolve any excursion beyond the sphere in O(1) draws with zero bias.

[[nodiscard]] inline double return_probability(double distance, double radius) {
    return radius / distance;
}

// The cosine of the landing point's angle against the axis from the sphere's
// center toward the walker, from one uniform unit draw. The draw fixes the
// straight-line separation between the walker and its landing point; the
// cosine follows from the triangle walker-center-landing.
[[nodiscard]] inline double harmonic_hit_cosine(double distance, double radius, double unit_draw) {
    double landing_separation =
        (distance * distance - radius * radius) /
        (2.0 * radius * unit_draw + distance - radius);

    return (
        distance * distance + radius * radius -
        landing_separation * landing_separation
    ) / (2.0 * distance * radius);
}

} // namespace geometry_art::dla

#endif //GEOMETRY_ART_DLA_SPAWN_SHELL_HARMONIC_HIT_HPP_
