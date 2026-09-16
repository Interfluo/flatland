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

# The C ABI test is written in C, so building it needs a C compiler that the
# rest of the project does not. Off by default; CI turns it on where the shell
# suite cannot build it itself (Windows, which has neither make nor a POSIX cc).
option(FLATLAND_BUILD_CAPI_TEST "Build and register the C ABI test program (needs a C compiler)" OFF)

# std::thread needs the platform threading library (e.g. -pthread on glibc).
find_package(Threads REQUIRED)

set(FLATLAND_CORE_SOURCES
    src/fl_mesh.cpp
    src/fl_deflate.cpp
    src/fl_parallel.cpp
    src/fl_io.cpp
    src/fl_field.cpp
    src/fl_raster.cpp
    src/fl_project.cpp
    src/fl_image.cpp
    src/fl_batch.cpp)

# ---------------------------------------------------------------------------
# Compiler settings that are about portability rather than taste. PRIVATE
# throughout: none of this is imposed on a consumer of the installed library.
# ---------------------------------------------------------------------------
function(flatland_target_flags tgt)
    if(MSVC)
        target_compile_options(${tgt} PRIVATE
            # The sources are UTF-8 (the comments carry the integral sign the
            # documentation uses). Without /utf-8 MSVC decodes them in the
            # machine's active code page and warns C4819 on every such file.
            /utf-8
            # The Makefile builds with -Wall; /W4 is the nearest MSVC setting.
            # Deliberately NOT /WX: the size_t/int narrowings MSVC flags here
            # are the same ones GCC's -Wconversion flags, and they are checked
            # and benign, so they should not be able to break a build.
            /W4)
        target_compile_definitions(${tgt} PRIVATE
            # windows.h, reachable through some toolchain headers, defines min
            # and max as macros, which breaks std::numeric_limits<T>::max().
            NOMINMAX
            WIN32_LEAN_AND_MEAN
            _CRT_SECURE_NO_WARNINGS)
    else()
        target_compile_options(${tgt} PRIVATE -Wall)
    endif()
endfunction()

# ---------------------------------------------------------------------------
# The engine, as an object library so the CLI and both libraries share one
# compile. Note that capi.cpp is NOT here: it is the only file whose object
# code differs between the two libraries (see FL_API below), so each library
# compiles its own copy.
# ---------------------------------------------------------------------------
add_library(flatland_core OBJECT ${FLATLAND_CORE_SOURCES})
target_include_directories(flatland_core PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(flatland_core PUBLIC Threads::Threads)
set_target_properties(flatland_core PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN ON)
flatland_target_flags(flatland_core)

# ---------------------------------------------------------------------------
# Libraries. Both carry the same C ABI; only flatland.h is public.
#
# On Windows, FL_API is __declspec(dllexport) in the DLL and nothing at all in
# the static library, so capi.cpp is compiled once per library with the matching
# macro. Compiling it once and archiving the same object into both — which is
# what an object library would do — stamps dllexport on the static library too,
# and then every program that links the .lib re-exports the whole C ABI.
#
# Exports are marked explicitly in the header, so WINDOWS_EXPORT_ALL_SYMBOLS is
# deliberately left off: turning it on would export the engine's internal C++
# symbols as well and tie consumers to this exact compiler and STL version.
# ---------------------------------------------------------------------------
add_library(flatland_shared SHARED src/capi.cpp)
target_compile_definitions(flatland_shared PRIVATE FLATLAND_BUILD_SHARED=1)

add_library(flatland_static STATIC src/capi.cpp)
# PUBLIC, not PRIVATE: a consumer of the static library has to define
# FLATLAND_STATIC as well, or flatland.h hands it dllimport declarations and the
# link fails looking for __imp_fl_*. Carrying it on the target's interface means
# anyone linking this from CMake gets it for free.
target_compile_definitions(flatland_static PUBLIC FLATLAND_STATIC=1)

foreach(tgt flatland_shared flatland_static)
    target_link_libraries(${tgt} PRIVATE flatland_core)
    target_link_libraries(${tgt} PUBLIC Threads::Threads)
    target_include_directories(${tgt} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
    target_include_directories(${tgt} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>)
    set_target_properties(${tgt} PROPERTIES
        OUTPUT_NAME flatland
        POSITION_INDEPENDENT_CODE ON
        CXX_VISIBILITY_PRESET hidden
        VISIBILITY_INLINES_HIDDEN ON)
    flatland_target_flags(${tgt})
endforeach()

if(MSVC)
    # The static library and the DLL's import library would both be called
    # flatland.lib and land in the same directory, so one silently overwrites
    # the other. Only MSVC has this collision: MinGW names the import library
    # libflatland.dll.a, which does not clash with libflatland.a.
    set_target_properties(flatland_static PROPERTIES OUTPUT_NAME flatland_static)
endif()

if(NOT WIN32)
    # PE has no soname. Setting VERSION on Windows only renames the DLL
    # (libflatland-3.dll under MinGW), which breaks the expected file name.
    set_target_properties(flatland_shared PROPERTIES
        VERSION   ${PROJECT_VERSION}
        SOVERSION ${PROJECT_VERSION_MAJOR})
endif()

add_library(FlatLand::flatland        ALIAS flatland_shared)
add_library(FlatLand::flatland_static ALIAS flatland_static)

# ---------------------------------------------------------------------------
# Command-line front end.
# ---------------------------------------------------------------------------
add_executable(flatland cli/main.cpp)
target_link_libraries(flatland PRIVATE flatland_core)
flatland_target_flags(flatland)

# ---------------------------------------------------------------------------
# Tests: the shell suite drives the CLI and, on POSIX, builds the C ABI test
# itself with make and cc. Windows has neither, so there the C ABI test is
# built here instead and the shell suite runs its other four suites.
# ---------------------------------------------------------------------------
include(CTest)
if(BUILD_TESTING)
    if(WIN32)
        # 'capi' shells out to make and cc and looks for libflatland.so;
        # FLATLAND_BUILD_CAPI_TEST covers the same ground natively.
        set(FLATLAND_SUITES geometry parsing cli batch)
    else()
        set(FLATLAND_SUITES "")     # empty == every suite
    endif()

    # Prefer Git for Windows' bash over the WSL stub in System32, which is on
    # PATH ahead of it and cannot see the build tree.
    find_program(FLATLAND_BASH NAMES bash
                 HINTS "$ENV{ProgramFiles}/Git/bin" "$ENV{ProgramFiles}/Git/usr/bin")
    if(FLATLAND_BASH)
        add_test(NAME suite
                 COMMAND ${FLATLAND_BASH} ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_tests.sh
                         $<TARGET_FILE:flatland> ${FLATLAND_SUITES}
                 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    else()
        message(STATUS "bash not found; the shell test suite will not be registered")
    endif()
endif()

if(FLATLAND_BUILD_CAPI_TEST)
    enable_language(C)

    # One field value per cube vertex; the values themselves are arbitrary, the
    # row count is what the test checks against the mesh.
    file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/examples/cube_area/cube.obj"
         FLATLAND_CUBE_VERTS REGEX "^v ")
    set(FLATLAND_CAPI_FIELD "")
    set(FLATLAND_CAPI_N 0)
    foreach(line IN LISTS FLATLAND_CUBE_VERTS)
        math(EXPR FLATLAND_CAPI_N "(${FLATLAND_CAPI_N} + 1) % 7")
        string(APPEND FLATLAND_CAPI_FIELD "${FLATLAND_CAPI_N}\n")
    endforeach()
    file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/capi_field.txt" "${FLATLAND_CAPI_FIELD}")

    # Both link modes are built: a missing export only shows up in the shared
    # case, and a stray dllexport only shows up in the static one.
    foreach(link shared static)
        add_executable(capi_test_${link} tests/test_capi.c)
        set_target_properties(capi_test_${link} PROPERTIES C_STANDARD 99)
        target_link_libraries(capi_test_${link} PRIVATE flatland_${link})
        if(BUILD_TESTING)
            add_test(NAME capi_${link}
                     COMMAND capi_test_${link}
                             ${CMAKE_CURRENT_SOURCE_DIR}/examples/cube_area/cube.obj
                             ${CMAKE_CURRENT_BINARY_DIR}/capi_field.txt)
        endif()
    endforeach()

    if(WIN32)
        # The DLL has to sit next to the executable that loads it; CMake puts
        # both in the same runtime output directory, but a multi-config
        # generator only guarantees that per configuration, so copy explicitly.
        add_custom_command(TARGET capi_test_shared POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                    $<TARGET_FILE:flatland_shared> $<TARGET_FILE_DIR:capi_test_shared>)
    elseif(NOT APPLE)
        set_target_properties(capi_test_shared PROPERTIES
            BUILD_RPATH $<TARGET_FILE_DIR:flatland_shared>)
    endif()
endif()

install(TARGETS flatland flatland_shared flatland_static
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
install(FILES include/flatland.h DESTINATION include)
