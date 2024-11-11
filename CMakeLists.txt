# Minimum version of CMake required
cmake_minimum_required(VERSION 3.10)

# Project name and version
project(pimavilo_mapper VERSION 1.0.0)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#Add subdirectory bioparser
add_subdirectory(external/bioparser)

# Add the executable
add_executable(pimavilo_mapper src/pimavilo_mapper.cpp)

#Link bioparser to the executable
target_link_libraries(pimavilo_mapper PRIVATE bioparser)

# Add version information
target_compile_definitions(pimavilo_mapper PRIVATE PROJECT_VERSION="${PROJECT_VERSION}")
