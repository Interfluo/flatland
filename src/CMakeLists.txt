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

add_executable(flatland main.cpp)

# std::thread needs the platform threading library (e.g. -pthread on glibc).
find_package(Threads REQUIRED)
target_link_libraries(flatland PRIVATE Threads::Threads)

install(TARGETS flatland DESTINATION bin)
