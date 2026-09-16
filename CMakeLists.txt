cmake_minimum_required(VERSION 3.14)

project(FlatLand
    VERSION 2.0
    DESCRIPTION "Projected-area and field analysis of arbitrary meshes"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)        # strict ISO C++ for portability

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Library half: everything except the command-line front end.
add_library(flatland_core STATIC
    src/fl_mesh.cpp
    src/fl_io.cpp
    src/fl_field.cpp
    src/fl_raster.cpp
    src/fl_project.cpp
    src/fl_image.cpp
    src/fl_batch.cpp)
target_include_directories(flatland_core PUBLIC src)

add_executable(flatland cli/main.cpp)
target_link_libraries(flatland PRIVATE flatland_core)

# std::thread needs the platform threading library (e.g. -pthread on glibc).
find_package(Threads REQUIRED)
target_link_libraries(flatland_core PUBLIC Threads::Threads)

install(TARGETS flatland DESTINATION bin)
