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
include test/CMakeFiles/test_fft.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include test/CMakeFiles/test_fft.dir/compiler_depend.make

# Include the progress variables for this target.
include test/CMakeFiles/test_fft.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/test_fft.dir/flags.make

test/CMakeFiles/test_fft.dir/test_fft.cpp.o: test/CMakeFiles/test_fft.dir/flags.make
test/CMakeFiles/test_fft.dir/test_fft.cpp.o: test/test_fft.cpp
test/CMakeFiles/test_fft.dir/test_fft.cpp.o: test/CMakeFiles/test_fft.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jiamian/rtfmm/3drpp_appx/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/CMakeFiles/test_fft.dir/test_fft.cpp.o"
	cd /home/jiamian/rtfmm/3drpp_appx/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT test/CMakeFiles/test_fft.dir/test_fft.cpp.o -MF CMakeFiles/test_fft.dir/test_fft.cpp.o.d -o CMakeFiles/test_fft.dir/test_fft.cpp.o -c /home/jiamian/rtfmm/3drpp_appx/test/test_fft.cpp

test/CMakeFiles/test_fft.dir/test_fft.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_fft.dir/test_fft.cpp.i"
	cd /home/jiamian/rtfmm/3drpp_appx/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jiamian/rtfmm/3drpp_appx/test/test_fft.cpp > CMakeFiles/test_fft.dir/test_fft.cpp.i

test/CMakeFiles/test_fft.dir/test_fft.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_fft.dir/test_fft.cpp.s"
	cd /home/jiamian/rtfmm/3drpp_appx/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jiamian/rtfmm/3drpp_appx/test/test_fft.cpp -o CMakeFiles/test_fft.dir/test_fft.cpp.s

# Object files for target test_fft
test_fft_OBJECTS = \
"CMakeFiles/test_fft.dir/test_fft.cpp.o"

# External object files for target test_fft
test_fft_EXTERNAL_OBJECTS =

test/test_fft: test/CMakeFiles/test_fft.dir/test_fft.cpp.o
test/test_fft: test/CMakeFiles/test_fft.dir/build.make
test/test_fft: src/libcore.a
test/test_fft: /usr/lib/gcc/x86_64-linux-gnu/9/libgomp.so
test/test_fft: /usr/lib/x86_64-linux-gnu/libpthread.so
test/test_fft: /usr/lib/x86_64-linux-gnu/libopenblas.so
test/test_fft: test/CMakeFiles/test_fft.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/jiamian/rtfmm/3drpp_appx/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable test_fft"
	cd /home/jiamian/rtfmm/3drpp_appx/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_fft.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/test_fft.dir/build: test/test_fft
.PHONY : test/CMakeFiles/test_fft.dir/build

test/CMakeFiles/test_fft.dir/clean:
	cd /home/jiamian/rtfmm/3drpp_appx/test && $(CMAKE_COMMAND) -P CMakeFiles/test_fft.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/test_fft.dir/clean

test/CMakeFiles/test_fft.dir/depend:
	cd /home/jiamian/rtfmm/3drpp_appx && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jiamian/rtfmm/3drpp_appx /home/jiamian/rtfmm/3drpp_appx/test /home/jiamian/rtfmm/3drpp_appx /home/jiamian/rtfmm/3drpp_appx/test /home/jiamian/rtfmm/3drpp_appx/test/CMakeFiles/test_fft.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/test_fft.dir/depend

