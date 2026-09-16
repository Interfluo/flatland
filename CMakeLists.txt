cmake_minimum_required(VERSION 3.14)

project(FlatLand
    VERSION 3.0.0
    DESCRIPTION "Projected area and scalar-field integrals over arbitrary meshes"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)        # strict ISO C++ for portability

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release)
endif()

# std::thread needs the platform threading library (e.g. -pthread on glibc).
find_package(Threads REQUIRED)

set(FLATLAND_CORE_SOURCES
    src/fl_mesh.cpp
    src/fl_parallel.cpp
    src/fl_io.cpp
    src/fl_field.cpp
    src/fl_raster.cpp
    src/fl_project.cpp
    src/fl_image.cpp
    src/fl_batch.cpp)

# ---------------------------------------------------------------------------
# The engine, as an object library so the static and shared libraries and the
# CLI can all share one compile without duplicating the source list.
# ---------------------------------------------------------------------------
add_library(flatland_objects OBJECT ${FLATLAND_CORE_SOURCES} src/capi.cpp)
target_include_directories(flatland_objects PUBLIC src include)
target_link_libraries(flatland_objects PUBLIC Threads::Threads)
set_target_properties(flatland_objects PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN ON)
target_compile_definitions(flatland_objects PRIVATE FLATLAND_BUILD_SHARED=1)

# ---------------------------------------------------------------------------
# Libraries. Both carry the same C ABI; only flatland.h is public.
# ---------------------------------------------------------------------------
add_library(flatland_shared SHARED $<TARGET_OBJECTS:flatland_objects>)
add_library(flatland_static STATIC $<TARGET_OBJECTS:flatland_objects>)

foreach(tgt flatland_shared flatland_static)
    target_include_directories(${tgt} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>)
    target_link_libraries(${tgt} PUBLIC Threads::Threads)
    set_target_properties(${tgt} PROPERTIES OUTPUT_NAME flatland)
endforeach()

set_target_properties(flatland_shared PROPERTIES
    VERSION   ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR})

add_library(FlatLand::flatland ALIAS flatland_shared)

# ---------------------------------------------------------------------------
# Command-line front end.
# ---------------------------------------------------------------------------
add_executable(flatland cli/main.cpp ${FLATLAND_CORE_SOURCES})
target_include_directories(flatland PRIVATE src include)
target_link_libraries(flatland PRIVATE Threads::Threads)

# ---------------------------------------------------------------------------
# Tests: the shell suite drives the CLI and builds the C ABI test itself.
# ---------------------------------------------------------------------------
include(CTest)
if(BUILD_TESTING)
    add_test(NAME suite
             COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_tests.sh $<TARGET_FILE:flatland>
             WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif()

install(TARGETS flatland flatland_shared flatland_static
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
install(FILES include/flatland.h DESTINATION include)
