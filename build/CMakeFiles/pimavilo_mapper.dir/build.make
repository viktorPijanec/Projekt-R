# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /mnt/c/Users/Lovre/Git/Projekt-R

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /mnt/c/Users/Lovre/Git/Projekt-R/build

# Include any dependencies generated for this target.
include CMakeFiles/pimavilo_mapper.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/pimavilo_mapper.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/pimavilo_mapper.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/pimavilo_mapper.dir/flags.make

CMakeFiles/pimavilo_mapper.dir/src/pimavilo_mapper.cpp.o: CMakeFiles/pimavilo_mapper.dir/flags.make
CMakeFiles/pimavilo_mapper.dir/src/pimavilo_mapper.cpp.o: ../src/pimavilo_mapper.cpp
CMakeFiles/pimavilo_mapper.dir/src/pimavilo_mapper.cpp.o: CMakeFiles/pimavilo_mapper.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/c/Users/Lovre/Git/Projekt-R/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/pimavilo_mapper.dir/src/pimavilo_mapper.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/pimavilo_mapper.dir/src/pimavilo_mapper.cpp.o -MF CMakeFiles/pimavilo_mapper.dir/src/pimavilo_mapper.cpp.o.d -o CMakeFiles/pimavilo_mapper.dir/src/pimavilo_mapper.cpp.o -c /mnt/c/Users/Lovre/Git/Projekt-R/src/pimavilo_mapper.cpp

CMakeFiles/pimavilo_mapper.dir/src/pimavilo_mapper.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/pimavilo_mapper.dir/src/pimavilo_mapper.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/c/Users/Lovre/Git/Projekt-R/src/pimavilo_mapper.cpp > CMakeFiles/pimavilo_mapper.dir/src/pimavilo_mapper.cpp.i

CMakeFiles/pimavilo_mapper.dir/src/pimavilo_mapper.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/pimavilo_mapper.dir/src/pimavilo_mapper.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/c/Users/Lovre/Git/Projekt-R/src/pimavilo_mapper.cpp -o CMakeFiles/pimavilo_mapper.dir/src/pimavilo_mapper.cpp.s

# Object files for target pimavilo_mapper
pimavilo_mapper_OBJECTS = \
"CMakeFiles/pimavilo_mapper.dir/src/pimavilo_mapper.cpp.o"

# External object files for target pimavilo_mapper
pimavilo_mapper_EXTERNAL_OBJECTS =

pimavilo_mapper: CMakeFiles/pimavilo_mapper.dir/src/pimavilo_mapper.cpp.o
pimavilo_mapper: CMakeFiles/pimavilo_mapper.dir/build.make
pimavilo_mapper: /usr/lib/x86_64-linux-gnu/libz.so
pimavilo_mapper: CMakeFiles/pimavilo_mapper.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/mnt/c/Users/Lovre/Git/Projekt-R/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable pimavilo_mapper"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/pimavilo_mapper.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/pimavilo_mapper.dir/build: pimavilo_mapper
.PHONY : CMakeFiles/pimavilo_mapper.dir/build

CMakeFiles/pimavilo_mapper.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/pimavilo_mapper.dir/cmake_clean.cmake
.PHONY : CMakeFiles/pimavilo_mapper.dir/clean

CMakeFiles/pimavilo_mapper.dir/depend:
	cd /mnt/c/Users/Lovre/Git/Projekt-R/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /mnt/c/Users/Lovre/Git/Projekt-R /mnt/c/Users/Lovre/Git/Projekt-R /mnt/c/Users/Lovre/Git/Projekt-R/build /mnt/c/Users/Lovre/Git/Projekt-R/build /mnt/c/Users/Lovre/Git/Projekt-R/build/CMakeFiles/pimavilo_mapper.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/pimavilo_mapper.dir/depend

