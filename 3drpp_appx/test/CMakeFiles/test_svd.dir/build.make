# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.20

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
CMAKE_COMMAND = /mnt/nfs/packages/x86_64/cmake/cmake-3.20.0/bin/cmake

# The command to remove a file.
RM = /mnt/nfs/packages/x86_64/cmake/cmake-3.20.0/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/jiamian/rtfmm/3drpp_appx

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/jiamian/rtfmm/3drpp_appx

# Include any dependencies generated for this target.
include test/CMakeFiles/test_svd.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include test/CMakeFiles/test_svd.dir/compiler_depend.make

# Include the progress variables for this target.
include test/CMakeFiles/test_svd.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/test_svd.dir/flags.make

test/CMakeFiles/test_svd.dir/test_svd.cpp.o: test/CMakeFiles/test_svd.dir/flags.make
test/CMakeFiles/test_svd.dir/test_svd.cpp.o: test/test_svd.cpp
test/CMakeFiles/test_svd.dir/test_svd.cpp.o: test/CMakeFiles/test_svd.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jiamian/rtfmm/3drpp_appx/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/CMakeFiles/test_svd.dir/test_svd.cpp.o"
	cd /home/jiamian/rtfmm/3drpp_appx/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT test/CMakeFiles/test_svd.dir/test_svd.cpp.o -MF CMakeFiles/test_svd.dir/test_svd.cpp.o.d -o CMakeFiles/test_svd.dir/test_svd.cpp.o -c /home/jiamian/rtfmm/3drpp_appx/test/test_svd.cpp

test/CMakeFiles/test_svd.dir/test_svd.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_svd.dir/test_svd.cpp.i"
	cd /home/jiamian/rtfmm/3drpp_appx/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jiamian/rtfmm/3drpp_appx/test/test_svd.cpp > CMakeFiles/test_svd.dir/test_svd.cpp.i

test/CMakeFiles/test_svd.dir/test_svd.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_svd.dir/test_svd.cpp.s"
	cd /home/jiamian/rtfmm/3drpp_appx/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jiamian/rtfmm/3drpp_appx/test/test_svd.cpp -o CMakeFiles/test_svd.dir/test_svd.cpp.s

# Object files for target test_svd
test_svd_OBJECTS = \
"CMakeFiles/test_svd.dir/test_svd.cpp.o"

# External object files for target test_svd
test_svd_EXTERNAL_OBJECTS =

test/test_svd: test/CMakeFiles/test_svd.dir/test_svd.cpp.o
test/test_svd: test/CMakeFiles/test_svd.dir/build.make
test/test_svd: src/libcore.a
test/test_svd: /usr/lib/gcc/x86_64-linux-gnu/9/libgomp.so
test/test_svd: /usr/lib/x86_64-linux-gnu/libpthread.so
test/test_svd: /usr/lib/x86_64-linux-gnu/libopenblas.so
test/test_svd: test/CMakeFiles/test_svd.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/jiamian/rtfmm/3drpp_appx/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable test_svd"
	cd /home/jiamian/rtfmm/3drpp_appx/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_svd.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/test_svd.dir/build: test/test_svd
.PHONY : test/CMakeFiles/test_svd.dir/build

test/CMakeFiles/test_svd.dir/clean:
	cd /home/jiamian/rtfmm/3drpp_appx/test && $(CMAKE_COMMAND) -P CMakeFiles/test_svd.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/test_svd.dir/clean

test/CMakeFiles/test_svd.dir/depend:
	cd /home/jiamian/rtfmm/3drpp_appx && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jiamian/rtfmm/3drpp_appx /home/jiamian/rtfmm/3drpp_appx/test /home/jiamian/rtfmm/3drpp_appx /home/jiamian/rtfmm/3drpp_appx/test /home/jiamian/rtfmm/3drpp_appx/test/CMakeFiles/test_svd.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/test_svd.dir/depend

