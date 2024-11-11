# Specify the minimum version of CMake required
cmake_minimum_required(VERSION 3.10)

# Set the project name and version
project(pimavilo_mapper VERSION 1.0.0)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add the executable
add_executable(pimavilo_mapper src/pimavilo_mapper.cpp)

# Add version information
target_compile_definitions(pimavilo_mapper PRIVATE PROJECT_VERSION="${PROJECT_VERSION}")
