# Minimum version of CMake required
cmake_minimum_required(VERSION 3.10)

# Project name and version
project(pimavilo_mapper VERSION 1.0.0)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#Enable OpenMP
find_package(OpenMP REQUIRED)
if (OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif ()

#Add subdirectory bioparser
add_subdirectory(external/bioparser)

# Add the pimavilo library to pimavilo_alignment and pimavilo_minimizers
add_library(pimavilo_alignment STATIC
    src/pimavilo/Align.cpp
)
add_library(pimavilo_minimizers STATIC
    src/pimavilo/Minimizers.cpp
)

# Include the pimavilo library headers
target_include_directories(pimavilo_alignment
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/pimavilo
)
target_include_directories(pimavilo_minimizers
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/pimavilo
)

# Add the executable
add_executable(pimavilo_mapper src/pimavilo_mapper.cpp)

# Link bioparser and pimavilo to the executable
target_link_libraries(pimavilo_mapper
    PRIVATE bioparser
    PRIVATE pimavilo_alignment
    PRIVATE pimavilo_minimizers
    PRIVATE OpenMP::OpenMP_CXX
)
# Add version information
target_compile_definitions(pimavilo_mapper PRIVATE PROJECT_VERSION="${PROJECT_VERSION}")
