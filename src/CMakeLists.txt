cmake_minimum_required(VERSION 3.14)  # A reasonably modern version; works everywhere

project(FlatLand
    VERSION 0.1
    DESCRIPTION "Projected Areas of arbitrary meshesl"
    LANGUAGES CXX)

# Set the C++ standard you need (adjust if you're using features from newer standards)
set(CMAKE_CXX_STANDARD 17)          # Or 11, 14, 20, 23 depending on your code
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)       # Use strict ISO C++ (good for portability)

# Create the executable
add_executable(flatland main.cpp)     # "mytool" is the output binary name; change if you want

# Optional: Enable optimizations for release builds (default is fine otherwise)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Optional: Install target (copies binary to a standard location when you run `cmake --install`)
install(TARGETS flatland DESTINATION bin)
