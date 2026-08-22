cmake_minimum_required(VERSION 3.16...3.31)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(GlobeArt)

set(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE true)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_DISABLE_FIND_PACKAGE_WrapVulkanHeaders ON)

# Libraries

find_package(PkgConfig REQUIRED)
find_package(Qt6 REQUIRED COMPONENTS Core Gui OpenGLWidgets Svg)
find_package(CGAL REQUIRED COMPONENTS Qt6)
find_package(Boost REQUIRED CONFIG)
find_package(anl REQUIRED)
find_package(CLI11 REQUIRED CONFIG)
find_package(Eigen3 REQUIRED)
find_package(Threads REQUIRED)
find_path(LBFGSPP_INCLUDE_DIR LBFGS.h REQUIRED)

# Header-only core library

add_library(globe INTERFACE)
target_include_directories(globe INTERFACE ${CMAKE_SOURCE_DIR}/src ${LBFGSPP_INCLUDE_DIR})
target_link_libraries(globe INTERFACE CGAL::CGAL Eigen3::Eigen Boost::boost Threads::Threads)

add_library(globe_noise STATIC src/globe/fields/scalar/anl_implementation.cpp)
target_link_libraries(globe_noise PUBLIC globe anl)

# Executables

add_executable(generate_globe src/generate_globe.cpp)
target_link_libraries(generate_globe PRIVATE globe globe_noise CGAL::CGAL_Qt6 CLI11::CLI11)
target_compile_definitions(generate_globe PRIVATE CGAL_USE_BASIC_VIEWER CGAL_USE_BASIC_VIEWER_QT QT_NO_KEYWORDS)

add_executable(visualize_field src/visualize_field.cpp)
target_link_libraries(visualize_field PRIVATE globe CGAL::CGAL_Qt6)
target_compile_definitions(visualize_field PRIVATE CGAL_USE_BASIC_VIEWER CGAL_USE_BASIC_VIEWER_QT QT_NO_KEYWORDS)

add_executable(sphere_to_mesh src/sphere_to_mesh.cpp)
target_link_libraries(sphere_to_mesh PRIVATE globe CLI11::CLI11)

# Tests

pkg_search_module(GTEST REQUIRED gtest_main)
pkg_search_module(GMOCK REQUIRED gmock)

file(GLOB_RECURSE TEST_SOURCES src/globe/*_test.cpp test/*_test.cpp)

foreach (TEST_SRC ${TEST_SOURCES})
    file(RELATIVE_PATH TEST_REL_PATH ${CMAKE_SOURCE_DIR} ${TEST_SRC})
    string(REPLACE "src/globe/" "" TEST_REL_PATH ${TEST_REL_PATH})
    string(REPLACE "test/" "" TEST_REL_PATH ${TEST_REL_PATH})
    string(REPLACE "/" "_" TEST_NAME ${TEST_REL_PATH})
    string(REPLACE ".cpp" "" TEST_NAME ${TEST_NAME})
    add_executable(${TEST_NAME} ${TEST_SRC})
    target_link_libraries(${TEST_NAME} PRIVATE globe globe_noise ${GTEST_LDFLAGS} ${GMOCK_LDFLAGS})
    target_compile_options(${TEST_NAME} PRIVATE ${GTEST_CFLAGS} ${GMOCK_CFLAGS})
    target_compile_definitions(${TEST_NAME} PRIVATE
        GLOBE_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
        GLOBE_BINARY_DIR="${CMAKE_BINARY_DIR}")
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
endforeach ()

enable_testing()
