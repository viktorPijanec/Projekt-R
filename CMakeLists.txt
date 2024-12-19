# Minimum version of CMake required
cmake_minimum_required(VERSION 3.10)

# Project name and version
project(pimavilo_mapper VERSION 1.0.0)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add subdirectory for bioparser
add_subdirectory(external/bioparser)

# Add the pimavilo alignment library
add_library(pimavilo_alignment STATIC
    src/pimavilo/Align.cpp
)

# Add the pimavilo minimizers library
add_library(pimavilo_minimizers STATIC
    src/pimavilo/Minimizers.cpp
)

# Include directories for libraries
target_include_directories(pimavilo_alignment
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_include_directories(pimavilo_minimizers
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# Add the executable
add_executable(pimavilo_mapper src/pimavilo_mapper.cpp)

# Link bioparser, alignment, and minimizers libraries to the executable
target_link_libraries(pimavilo_mapper
    PRIVATE bioparser
    PRIVATE pimavilo_alignment
    PRIVATE pimavilo_minimizers
)

# Add version information
target_compile_definitions(pimavilo_mapper PRIVATE PROJECT_VERSION="${PROJECT_VERSION}")