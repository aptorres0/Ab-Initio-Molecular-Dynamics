# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.25

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
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.25.3/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.25.3/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/Users/ajpt/Documents/School/PHYS 7411 - Computational/dev/Project-3"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/Users/ajpt/Documents/School/PHYS 7411 - Computational/dev/Project-3/build"

# Include any dependencies generated for this target.
include CMakeFiles/main.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/main.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/main.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/main.dir/flags.make

CMakeFiles/main.dir/scripts/main.cxx.o: CMakeFiles/main.dir/flags.make
CMakeFiles/main.dir/scripts/main.cxx.o: /Users/ajpt/Documents/School/PHYS\ 7411\ -\ Computational/dev/Project-3/scripts/main.cxx
CMakeFiles/main.dir/scripts/main.cxx.o: CMakeFiles/main.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/Users/ajpt/Documents/School/PHYS 7411 - Computational/dev/Project-3/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/main.dir/scripts/main.cxx.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/main.dir/scripts/main.cxx.o -MF CMakeFiles/main.dir/scripts/main.cxx.o.d -o CMakeFiles/main.dir/scripts/main.cxx.o -c "/Users/ajpt/Documents/School/PHYS 7411 - Computational/dev/Project-3/scripts/main.cxx"

CMakeFiles/main.dir/scripts/main.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main.dir/scripts/main.cxx.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/Users/ajpt/Documents/School/PHYS 7411 - Computational/dev/Project-3/scripts/main.cxx" > CMakeFiles/main.dir/scripts/main.cxx.i

CMakeFiles/main.dir/scripts/main.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main.dir/scripts/main.cxx.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/Users/ajpt/Documents/School/PHYS 7411 - Computational/dev/Project-3/scripts/main.cxx" -o CMakeFiles/main.dir/scripts/main.cxx.s

# Object files for target main
main_OBJECTS = \
"CMakeFiles/main.dir/scripts/main.cxx.o"

# External object files for target main
main_EXTERNAL_OBJECTS =

scripts/main: CMakeFiles/main.dir/scripts/main.cxx.o
scripts/main: CMakeFiles/main.dir/build.make
scripts/main: libproject3.dylib
scripts/main: CMakeFiles/main.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/Users/ajpt/Documents/School/PHYS 7411 - Computational/dev/Project-3/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable scripts/main"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/main.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/main.dir/build: scripts/main
.PHONY : CMakeFiles/main.dir/build

CMakeFiles/main.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/main.dir/cmake_clean.cmake
.PHONY : CMakeFiles/main.dir/clean

CMakeFiles/main.dir/depend:
	cd "/Users/ajpt/Documents/School/PHYS 7411 - Computational/dev/Project-3/build" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/Users/ajpt/Documents/School/PHYS 7411 - Computational/dev/Project-3" "/Users/ajpt/Documents/School/PHYS 7411 - Computational/dev/Project-3" "/Users/ajpt/Documents/School/PHYS 7411 - Computational/dev/Project-3/build" "/Users/ajpt/Documents/School/PHYS 7411 - Computational/dev/Project-3/build" "/Users/ajpt/Documents/School/PHYS 7411 - Computational/dev/Project-3/build/CMakeFiles/main.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/main.dir/depend

