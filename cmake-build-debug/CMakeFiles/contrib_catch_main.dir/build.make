# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.17

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
CMAKE_COMMAND = /home/user/.local/share/JetBrains/Toolbox/apps/CLion/ch-0/202.7319.72/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /home/user/.local/share/JetBrains/Toolbox/apps/CLion/ch-0/202.7319.72/bin/cmake/linux/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/user/Mine/raytracer

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/user/Mine/raytracer/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/contrib_catch_main.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/contrib_catch_main.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/contrib_catch_main.dir/flags.make

CMakeFiles/contrib_catch_main.dir/contrib/catch/catch_main.cpp.o: CMakeFiles/contrib_catch_main.dir/flags.make
CMakeFiles/contrib_catch_main.dir/contrib/catch/catch_main.cpp.o: ../contrib/catch/catch_main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/user/Mine/raytracer/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/contrib_catch_main.dir/contrib/catch/catch_main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/contrib_catch_main.dir/contrib/catch/catch_main.cpp.o -c /home/user/Mine/raytracer/contrib/catch/catch_main.cpp

CMakeFiles/contrib_catch_main.dir/contrib/catch/catch_main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/contrib_catch_main.dir/contrib/catch/catch_main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/user/Mine/raytracer/contrib/catch/catch_main.cpp > CMakeFiles/contrib_catch_main.dir/contrib/catch/catch_main.cpp.i

CMakeFiles/contrib_catch_main.dir/contrib/catch/catch_main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/contrib_catch_main.dir/contrib/catch/catch_main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/user/Mine/raytracer/contrib/catch/catch_main.cpp -o CMakeFiles/contrib_catch_main.dir/contrib/catch/catch_main.cpp.s

# Object files for target contrib_catch_main
contrib_catch_main_OBJECTS = \
"CMakeFiles/contrib_catch_main.dir/contrib/catch/catch_main.cpp.o"

# External object files for target contrib_catch_main
contrib_catch_main_EXTERNAL_OBJECTS =

libcontrib_catch_main.a: CMakeFiles/contrib_catch_main.dir/contrib/catch/catch_main.cpp.o
libcontrib_catch_main.a: CMakeFiles/contrib_catch_main.dir/build.make
libcontrib_catch_main.a: CMakeFiles/contrib_catch_main.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/user/Mine/raytracer/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libcontrib_catch_main.a"
	$(CMAKE_COMMAND) -P CMakeFiles/contrib_catch_main.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/contrib_catch_main.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/contrib_catch_main.dir/build: libcontrib_catch_main.a

.PHONY : CMakeFiles/contrib_catch_main.dir/build

CMakeFiles/contrib_catch_main.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/contrib_catch_main.dir/cmake_clean.cmake
.PHONY : CMakeFiles/contrib_catch_main.dir/clean

CMakeFiles/contrib_catch_main.dir/depend:
	cd /home/user/Mine/raytracer/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/user/Mine/raytracer /home/user/Mine/raytracer /home/user/Mine/raytracer/cmake-build-debug /home/user/Mine/raytracer/cmake-build-debug /home/user/Mine/raytracer/cmake-build-debug/CMakeFiles/contrib_catch_main.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/contrib_catch_main.dir/depend

