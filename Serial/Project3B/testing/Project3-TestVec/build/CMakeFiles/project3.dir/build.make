# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.26

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
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.26.1/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.26.1/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/ajpt/School/Project3-TestVec

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/ajpt/School/Project3-TestVec/build

# Include any dependencies generated for this target.
include CMakeFiles/project3.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/project3.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/project3.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/project3.dir/flags.make

CMakeFiles/project3.dir/src/Sim.cxx.o: CMakeFiles/project3.dir/flags.make
CMakeFiles/project3.dir/src/Sim.cxx.o: /Users/ajpt/School/Project3-TestVec/src/Sim.cxx
CMakeFiles/project3.dir/src/Sim.cxx.o: CMakeFiles/project3.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/ajpt/School/Project3-TestVec/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/project3.dir/src/Sim.cxx.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/project3.dir/src/Sim.cxx.o -MF CMakeFiles/project3.dir/src/Sim.cxx.o.d -o CMakeFiles/project3.dir/src/Sim.cxx.o -c /Users/ajpt/School/Project3-TestVec/src/Sim.cxx

CMakeFiles/project3.dir/src/Sim.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/project3.dir/src/Sim.cxx.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/ajpt/School/Project3-TestVec/src/Sim.cxx > CMakeFiles/project3.dir/src/Sim.cxx.i

CMakeFiles/project3.dir/src/Sim.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/project3.dir/src/Sim.cxx.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/ajpt/School/Project3-TestVec/src/Sim.cxx -o CMakeFiles/project3.dir/src/Sim.cxx.s

# Object files for target project3
project3_OBJECTS = \
"CMakeFiles/project3.dir/src/Sim.cxx.o"

# External object files for target project3
project3_EXTERNAL_OBJECTS =

lib/libproject3.dylib: CMakeFiles/project3.dir/src/Sim.cxx.o
lib/libproject3.dylib: CMakeFiles/project3.dir/build.make
lib/libproject3.dylib: CMakeFiles/project3.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/ajpt/School/Project3-TestVec/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX shared library lib/libproject3.dylib"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/project3.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/project3.dir/build: lib/libproject3.dylib
.PHONY : CMakeFiles/project3.dir/build

CMakeFiles/project3.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/project3.dir/cmake_clean.cmake
.PHONY : CMakeFiles/project3.dir/clean

CMakeFiles/project3.dir/depend:
	cd /Users/ajpt/School/Project3-TestVec/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/ajpt/School/Project3-TestVec /Users/ajpt/School/Project3-TestVec /Users/ajpt/School/Project3-TestVec/build /Users/ajpt/School/Project3-TestVec/build /Users/ajpt/School/Project3-TestVec/build/CMakeFiles/project3.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/project3.dir/depend

