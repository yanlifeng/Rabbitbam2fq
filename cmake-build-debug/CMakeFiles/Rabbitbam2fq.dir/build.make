# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake

# The command to remove a file.
RM = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/ylf9811/CLionProjects/Rabbitbam2fq

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/ylf9811/CLionProjects/Rabbitbam2fq/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/Rabbitbam2fq.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/Rabbitbam2fq.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/Rabbitbam2fq.dir/flags.make

CMakeFiles/Rabbitbam2fq.dir/test.cpp.o: CMakeFiles/Rabbitbam2fq.dir/flags.make
CMakeFiles/Rabbitbam2fq.dir/test.cpp.o: ../test.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/ylf9811/CLionProjects/Rabbitbam2fq/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/Rabbitbam2fq.dir/test.cpp.o"
	/usr/local/Cellar/gcc@8/8.4.0_1/bin/g++-8  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Rabbitbam2fq.dir/test.cpp.o -c /Users/ylf9811/CLionProjects/Rabbitbam2fq/test.cpp

CMakeFiles/Rabbitbam2fq.dir/test.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Rabbitbam2fq.dir/test.cpp.i"
	/usr/local/Cellar/gcc@8/8.4.0_1/bin/g++-8 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/ylf9811/CLionProjects/Rabbitbam2fq/test.cpp > CMakeFiles/Rabbitbam2fq.dir/test.cpp.i

CMakeFiles/Rabbitbam2fq.dir/test.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Rabbitbam2fq.dir/test.cpp.s"
	/usr/local/Cellar/gcc@8/8.4.0_1/bin/g++-8 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/ylf9811/CLionProjects/Rabbitbam2fq/test.cpp -o CMakeFiles/Rabbitbam2fq.dir/test.cpp.s

# Object files for target Rabbitbam2fq
Rabbitbam2fq_OBJECTS = \
"CMakeFiles/Rabbitbam2fq.dir/test.cpp.o"

# External object files for target Rabbitbam2fq
Rabbitbam2fq_EXTERNAL_OBJECTS =

Rabbitbam2fq: CMakeFiles/Rabbitbam2fq.dir/test.cpp.o
Rabbitbam2fq: CMakeFiles/Rabbitbam2fq.dir/build.make
Rabbitbam2fq: CMakeFiles/Rabbitbam2fq.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/ylf9811/CLionProjects/Rabbitbam2fq/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable Rabbitbam2fq"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Rabbitbam2fq.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/Rabbitbam2fq.dir/build: Rabbitbam2fq

.PHONY : CMakeFiles/Rabbitbam2fq.dir/build

CMakeFiles/Rabbitbam2fq.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/Rabbitbam2fq.dir/cmake_clean.cmake
.PHONY : CMakeFiles/Rabbitbam2fq.dir/clean

CMakeFiles/Rabbitbam2fq.dir/depend:
	cd /Users/ylf9811/CLionProjects/Rabbitbam2fq/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/ylf9811/CLionProjects/Rabbitbam2fq /Users/ylf9811/CLionProjects/Rabbitbam2fq /Users/ylf9811/CLionProjects/Rabbitbam2fq/cmake-build-debug /Users/ylf9811/CLionProjects/Rabbitbam2fq/cmake-build-debug /Users/ylf9811/CLionProjects/Rabbitbam2fq/cmake-build-debug/CMakeFiles/Rabbitbam2fq.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/Rabbitbam2fq.dir/depend

