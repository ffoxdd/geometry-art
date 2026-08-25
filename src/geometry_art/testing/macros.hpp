#ifndef GEOMETRY_ART_TESTING_MACROS_HPP_
#define GEOMETRY_ART_TESTING_MACROS_HPP_

#include <cstdlib>
#include <gtest/gtest.h>

#define REQUIRE_EXPENSIVE() \
    do { \
        const char* env_var = std::getenv("EXPENSIVE"); \
        if (!env_var || env_var[0] == '\0') { \
            GTEST_SKIP() << "Skipping expensive test (set EXPENSIVE=1 to run)"; \
        } \
    } while (0)

#endif //GEOMETRY_ART_TESTING_MACROS_HPP_
