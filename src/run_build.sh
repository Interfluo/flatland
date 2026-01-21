#!/bin/bash

# Start a timer for the entire script
echo "Starting build process..."
time (

    # Create a build directory (keeps source clean)
    echo "Creating build directory..."
    mkdir -p build 

    cd build

    # Configure/generate build files for your platform
    echo "Configuring project with CMake..."
    time cmake ..
    echo "CMake configuration finished."

    # Actually compile
    echo "Starting compilation..."
    time cmake --build . --config Release
    echo "Compilation finished."

)
echo "Total build process finished."
