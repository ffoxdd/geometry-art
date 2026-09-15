cmake_minimum_required(VERSION 3.16...3.31)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(GeometryArt)

set(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE true)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_DISABLE_FIND_PACKAGE_WrapVulkanHeaders ON)

# Libraries

find_package(PkgConfig REQUIRED)
find_package(CGAL REQUIRED)
find_package(Boost REQUIRED CONFIG)
find_package(anl REQUIRED)
find_package(CLI11 REQUIRED CONFIG)
find_package(Eigen3 REQUIRED)
find_package(Threads REQUIRED)
find_package(simdjson REQUIRED CONFIG)
find_path(LBFGSPP_INCLUDE_DIR LBFGS.h REQUIRED)

# Header-only core library

add_library(geometry_art INTERFACE)
target_include_directories(geometry_art INTERFACE ${CMAKE_SOURCE_DIR}/src ${LBFGSPP_INCLUDE_DIR})
target_link_libraries(geometry_art INTERFACE CGAL::CGAL Eigen3::Eigen Boost::boost Threads::Threads simdjson::simdjson)

add_library(geometry_art_noise STATIC src/geometry_art/fields/scalar/anl_implementation.cpp)
target_link_libraries(geometry_art_noise PUBLIC geometry_art anl)

# Executables

add_executable(tessellate src/tessellate.cpp)
target_link_libraries(tessellate PRIVATE geometry_art geometry_art_noise CLI11::CLI11)

add_executable(aggregate src/aggregate.cpp)
target_link_libraries(aggregate PRIVATE geometry_art CLI11::CLI11)

add_executable(skeletonize src/skeletonize.cpp)
target_link_libraries(skeletonize PRIVATE geometry_art CLI11::CLI11)

# Tests

pkg_search_module(GTEST REQUIRED gtest_main)
pkg_search_module(GMOCK REQUIRED gmock)

file(GLOB_RECURSE TEST_SOURCES src/geometry_art/*_test.cpp test/*_test.cpp)

foreach (TEST_SRC ${TEST_SOURCES})
    file(RELATIVE_PATH TEST_REL_PATH ${CMAKE_SOURCE_DIR} ${TEST_SRC})
    string(REPLACE "src/geometry_art/" "" TEST_REL_PATH ${TEST_REL_PATH})
    string(REPLACE "test/" "" TEST_REL_PATH ${TEST_REL_PATH})
    string(REPLACE "/" "_" TEST_NAME ${TEST_REL_PATH})
    string(REPLACE ".cpp" "" TEST_NAME ${TEST_NAME})
    add_executable(${TEST_NAME} ${TEST_SRC})
    target_link_libraries(${TEST_NAME} PRIVATE geometry_art geometry_art_noise ${GTEST_LDFLAGS} ${GMOCK_LDFLAGS})
    target_compile_options(${TEST_NAME} PRIVATE ${GTEST_CFLAGS} ${GMOCK_CFLAGS})
    target_compile_definitions(${TEST_NAME} PRIVATE
        GEOMETRY_ART_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
        GEOMETRY_ART_BINARY_DIR="${CMAKE_BINARY_DIR}")
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
endforeach ()

enable_testing()
