# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.14

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
CMAKE_COMMAND = /usr/bin/cmake3

# The command to remove a file.
RM = /usr/bin/cmake3 -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /u/re/alspellm/work/src/hpstr

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /u/re/alspellm/work/src/hpstr/build

# Include any dependencies generated for this target.
include processors/CMakeFiles/processors.dir/depend.make

# Include the progress variables for this target.
include processors/CMakeFiles/processors.dir/progress.make

# Include the compile flags for this target's objects.
include processors/CMakeFiles/processors.dir/flags.make

processors/CMakeFiles/processors.dir/src/ClusterOnTrackAnaProcessor.cxx.o: processors/CMakeFiles/processors.dir/flags.make
processors/CMakeFiles/processors.dir/src/ClusterOnTrackAnaProcessor.cxx.o: ../processors/src/ClusterOnTrackAnaProcessor.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/u/re/alspellm/work/src/hpstr/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object processors/CMakeFiles/processors.dir/src/ClusterOnTrackAnaProcessor.cxx.o"
	cd /u/re/alspellm/work/src/hpstr/build/processors && /opt/rh/devtoolset-8/root/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/processors.dir/src/ClusterOnTrackAnaProcessor.cxx.o -c /u/re/alspellm/work/src/hpstr/processors/src/ClusterOnTrackAnaProcessor.cxx

processors/CMakeFiles/processors.dir/src/ClusterOnTrackAnaProcessor.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/processors.dir/src/ClusterOnTrackAnaProcessor.cxx.i"
	cd /u/re/alspellm/work/src/hpstr/build/processors && /opt/rh/devtoolset-8/root/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /u/re/alspellm/work/src/hpstr/processors/src/ClusterOnTrackAnaProcessor.cxx > CMakeFiles/processors.dir/src/ClusterOnTrackAnaProcessor.cxx.i

processors/CMakeFiles/processors.dir/src/ClusterOnTrackAnaProcessor.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/processors.dir/src/ClusterOnTrackAnaProcessor.cxx.s"
	cd /u/re/alspellm/work/src/hpstr/build/processors && /opt/rh/devtoolset-8/root/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /u/re/alspellm/work/src/hpstr/processors/src/ClusterOnTrackAnaProcessor.cxx -o CMakeFiles/processors.dir/src/ClusterOnTrackAnaProcessor.cxx.s

processors/CMakeFiles/processors.dir/src/ECalDataProcessor.cxx.o: processors/CMakeFiles/processors.dir/flags.make
processors/CMakeFiles/processors.dir/src/ECalDataProcessor.cxx.o: ../processors/src/ECalDataProcessor.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/u/re/alspellm/work/src/hpstr/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object processors/CMakeFiles/processors.dir/src/ECalDataProcessor.cxx.o"
	cd /u/re/alspellm/work/src/hpstr/build/processors && /opt/rh/devtoolset-8/root/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/processors.dir/src/ECalDataProcessor.cxx.o -c /u/re/alspellm/work/src/hpstr/processors/src/ECalDataProcessor.cxx

processors/CMakeFiles/processors.dir/src/ECalDataProcessor.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/processors.dir/src/ECalDataProcessor.cxx.i"
	cd /u/re/alspellm/work/src/hpstr/build/processors && /opt/rh/devtoolset-8/root/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /u/re/alspellm/work/src/hpstr/processors/src/ECalDataProcessor.cxx > CMakeFiles/processors.dir/src/ECalDataProcessor.cxx.i

processors/CMakeFiles/processors.dir/src/ECalDataProcessor.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/processors.dir/src/ECalDataProcessor.cxx.s"
	cd /u/re/alspellm/work/src/hpstr/build/processors && /opt/rh/devtoolset-8/root/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /u/re/alspellm/work/src/hpstr/processors/src/ECalDataProcessor.cxx -o CMakeFiles/processors.dir/src/ECalDataProcessor.cxx.s

processors/CMakeFiles/processors.dir/src/EventProcessor.cxx.o: processors/CMakeFiles/processors.dir/flags.make
processors/CMakeFiles/processors.dir/src/EventProcessor.cxx.o: ../processors/src/EventProcessor.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/u/re/alspellm/work/src/hpstr/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object processors/CMakeFiles/processors.dir/src/EventProcessor.cxx.o"
	cd /u/re/alspellm/work/src/hpstr/build/processors && /opt/rh/devtoolset-8/root/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/processors.dir/src/EventProcessor.cxx.o -c /u/re/alspellm/work/src/hpstr/processors/src/EventProcessor.cxx

processors/CMakeFiles/processors.dir/src/EventProcessor.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/processors.dir/src/EventProcessor.cxx.i"
	cd /u/re/alspellm/work/src/hpstr/build/processors && /opt/rh/devtoolset-8/root/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /u/re/alspellm/work/src/hpstr/processors/src/EventProcessor.cxx > CMakeFiles/processors.dir/src/EventProcessor.cxx.i

processors/CMakeFiles/processors.dir/src/EventProcessor.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/processors.dir/src/EventProcessor.cxx.s"
	cd /u/re/alspellm/work/src/hpstr/build/processors && /opt/rh/devtoolset-8/root/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /u/re/alspellm/work/src/hpstr/processors/src/EventProcessor.cxx -o CMakeFiles/processors.dir/src/EventProcessor.cxx.s

processors/CMakeFiles/processors.dir/src/HPSEventProcessor.cxx.o: processors/CMakeFiles/processors.dir/flags.make
processors/CMakeFiles/processors.dir/src/HPSEventProcessor.cxx.o: ../processors/src/HPSEventProcessor.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/u/re/alspellm/work/src/hpstr/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object processors/CMakeFiles/processors.dir/src/HPSEventProcessor.cxx.o"
	cd /u/re/alspellm/work/src/hpstr/build/processors && /opt/rh/devtoolset-8/root/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/processors.dir/src/HPSEventProcessor.cxx.o -c /u/re/alspellm/work/src/hpstr/processors/src/HPSEventProcessor.cxx

processors/CMakeFiles/processors.dir/src/HPSEventProcessor.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/processors.dir/src/HPSEventProcessor.cxx.i"
	cd /u/re/alspellm/work/src/hpstr/build/processors && /opt/rh/devtoolset-8/root/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /u/re/alspellm/work/src/hpstr/processors/src/HPSEventProcessor.cxx > CMakeFiles/processors.dir/src/HPSEventProcessor.cxx.i

processors/CMakeFiles/processors.dir/src/HPSEventProcessor.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/processors.dir/src/HPSEventProcessor.cxx.s"
	cd /u/re/alspellm/work/src/hpstr/build/processors && /opt/rh/devtoolset-8/root/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /u/re/alspellm/work/src/hpstr/processors/src/HPSEventProcessor.cxx -o CMakeFiles/processors.dir/src/HPSEventProcessor.cxx.s

processors/CMakeFiles/processors.dir/src/ParticleProcessor.cxx.o: processors/CMakeFiles/processors.dir/flags.make
processors/CMakeFiles/processors.dir/src/ParticleProcessor.cxx.o: ../processors/src/ParticleProcessor.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/u/re/alspellm/work/src/hpstr/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object processors/CMakeFiles/processors.dir/src/ParticleProcessor.cxx.o"
	cd /u/re/alspellm/work/src/hpstr/build/processors && /opt/rh/devtoolset-8/root/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/processors.dir/src/ParticleProcessor.cxx.o -c /u/re/alspellm/work/src/hpstr/processors/src/ParticleProcessor.cxx

processors/CMakeFiles/processors.dir/src/ParticleProcessor.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/processors.dir/src/ParticleProcessor.cxx.i"
	cd /u/re/alspellm/work/src/hpstr/build/processors && /opt/rh/devtoolset-8/root/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /u/re/alspellm/work/src/hpstr/processors/src/ParticleProcessor.cxx > CMakeFiles/processors.dir/src/ParticleProcessor.cxx.i

processors/CMakeFiles/processors.dir/src/ParticleProcessor.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/processors.dir/src/ParticleProcessor.cxx.s"
	cd /u/re/alspellm/work/src/hpstr/build/processors && /opt/rh/devtoolset-8/root/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /u/re/alspellm/work/src/hpstr/processors/src/ParticleProcessor.cxx -o CMakeFiles/processors.dir/src/ParticleProcessor.cxx.s

processors/CMakeFiles/processors.dir/src/RefittedTracksProcessor.cxx.o: processors/CMakeFiles/processors.dir/flags.make
processors/CMakeFiles/processors.dir/src/RefittedTracksProcessor.cxx.o: ../processors/src/RefittedTracksProcessor.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/u/re/alspellm/work/src/hpstr/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object processors/CMakeFiles/processors.dir/src/RefittedTracksProcessor.cxx.o"
	cd /u/re/alspellm/work/src/hpstr/build/processors && /opt/rh/devtoolset-8/root/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/processors.dir/src/RefittedTracksProcessor.cxx.o -c /u/re/alspellm/work/src/hpstr/processors/src/RefittedTracksProcessor.cxx

processors/CMakeFiles/processors.dir/src/RefittedTracksProcessor.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/processors.dir/src/RefittedTracksProcessor.cxx.i"
	cd /u/re/alspellm/work/src/hpstr/build/processors && /opt/rh/devtoolset-8/root/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /u/re/alspellm/work/src/hpstr/processors/src/RefittedTracksProcessor.cxx > CMakeFiles/processors.dir/src/RefittedTracksProcessor.cxx.i

processors/CMakeFiles/processors.dir/src/RefittedTracksProcessor.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/processors.dir/src/RefittedTracksProcessor.cxx.s"
	cd /u/re/alspellm/work/src/hpstr/build/processors && /opt/rh/devtoolset-8/root/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /u/re/alspellm/work/src/hpstr/processors/src/RefittedTracksProcessor.cxx -o CMakeFiles/processors.dir/src/RefittedTracksProcessor.cxx.s

processors/CMakeFiles/processors.dir/src/SvtCondAnaProcessor.cxx.o: processors/CMakeFiles/processors.dir/flags.make
processors/CMakeFiles/processors.dir/src/SvtCondAnaProcessor.cxx.o: ../processors/src/SvtCondAnaProcessor.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/u/re/alspellm/work/src/hpstr/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object processors/CMakeFiles/processors.dir/src/SvtCondAnaProcessor.cxx.o"
	cd /u/re/alspellm/work/src/hpstr/build/processors && /opt/rh/devtoolset-8/root/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/processors.dir/src/SvtCondAnaProcessor.cxx.o -c /u/re/alspellm/work/src/hpstr/processors/src/SvtCondAnaProcessor.cxx

processors/CMakeFiles/processors.dir/src/SvtCondAnaProcessor.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/processors.dir/src/SvtCondAnaProcessor.cxx.i"
	cd /u/re/alspellm/work/src/hpstr/build/processors && /opt/rh/devtoolset-8/root/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /u/re/alspellm/work/src/hpstr/processors/src/SvtCondAnaProcessor.cxx > CMakeFiles/processors.dir/src/SvtCondAnaProcessor.cxx.i

processors/CMakeFiles/processors.dir/src/SvtCondAnaProcessor.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/processors.dir/src/SvtCondAnaProcessor.cxx.s"
	cd /u/re/alspellm/work/src/hpstr/build/processors && /opt/rh/devtoolset-8/root/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /u/re/alspellm/work/src/hpstr/processors/src/SvtCondAnaProcessor.cxx -o CMakeFiles/processors.dir/src/SvtCondAnaProcessor.cxx.s

processors/CMakeFiles/processors.dir/src/SvtDataProcessor.cxx.o: processors/CMakeFiles/processors.dir/flags.make
processors/CMakeFiles/processors.dir/src/SvtDataProcessor.cxx.o: ../processors/src/SvtDataProcessor.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/u/re/alspellm/work/src/hpstr/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object processors/CMakeFiles/processors.dir/src/SvtDataProcessor.cxx.o"
	cd /u/re/alspellm/work/src/hpstr/build/processors && /opt/rh/devtoolset-8/root/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/processors.dir/src/SvtDataProcessor.cxx.o -c /u/re/alspellm/work/src/hpstr/processors/src/SvtDataProcessor.cxx

processors/CMakeFiles/processors.dir/src/SvtDataProcessor.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/processors.dir/src/SvtDataProcessor.cxx.i"
	cd /u/re/alspellm/work/src/hpstr/build/processors && /opt/rh/devtoolset-8/root/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /u/re/alspellm/work/src/hpstr/processors/src/SvtDataProcessor.cxx > CMakeFiles/processors.dir/src/SvtDataProcessor.cxx.i

processors/CMakeFiles/processors.dir/src/SvtDataProcessor.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/processors.dir/src/SvtDataProcessor.cxx.s"
	cd /u/re/alspellm/work/src/hpstr/build/processors && /opt/rh/devtoolset-8/root/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /u/re/alspellm/work/src/hpstr/processors/src/SvtDataProcessor.cxx -o CMakeFiles/processors.dir/src/SvtDataProcessor.cxx.s

processors/CMakeFiles/processors.dir/src/SvtRawDataProcessor.cxx.o: processors/CMakeFiles/processors.dir/flags.make
processors/CMakeFiles/processors.dir/src/SvtRawDataProcessor.cxx.o: ../processors/src/SvtRawDataProcessor.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/u/re/alspellm/work/src/hpstr/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object processors/CMakeFiles/processors.dir/src/SvtRawDataProcessor.cxx.o"
	cd /u/re/alspellm/work/src/hpstr/build/processors && /opt/rh/devtoolset-8/root/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/processors.dir/src/SvtRawDataProcessor.cxx.o -c /u/re/alspellm/work/src/hpstr/processors/src/SvtRawDataProcessor.cxx

processors/CMakeFiles/processors.dir/src/SvtRawDataProcessor.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/processors.dir/src/SvtRawDataProcessor.cxx.i"
	cd /u/re/alspellm/work/src/hpstr/build/processors && /opt/rh/devtoolset-8/root/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /u/re/alspellm/work/src/hpstr/processors/src/SvtRawDataProcessor.cxx > CMakeFiles/processors.dir/src/SvtRawDataProcessor.cxx.i

processors/CMakeFiles/processors.dir/src/SvtRawDataProcessor.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/processors.dir/src/SvtRawDataProcessor.cxx.s"
	cd /u/re/alspellm/work/src/hpstr/build/processors && /opt/rh/devtoolset-8/root/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /u/re/alspellm/work/src/hpstr/processors/src/SvtRawDataProcessor.cxx -o CMakeFiles/processors.dir/src/SvtRawDataProcessor.cxx.s

processors/CMakeFiles/processors.dir/src/TrackingProcessor.cxx.o: processors/CMakeFiles/processors.dir/flags.make
processors/CMakeFiles/processors.dir/src/TrackingProcessor.cxx.o: ../processors/src/TrackingProcessor.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/u/re/alspellm/work/src/hpstr/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object processors/CMakeFiles/processors.dir/src/TrackingProcessor.cxx.o"
	cd /u/re/alspellm/work/src/hpstr/build/processors && /opt/rh/devtoolset-8/root/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/processors.dir/src/TrackingProcessor.cxx.o -c /u/re/alspellm/work/src/hpstr/processors/src/TrackingProcessor.cxx

processors/CMakeFiles/processors.dir/src/TrackingProcessor.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/processors.dir/src/TrackingProcessor.cxx.i"
	cd /u/re/alspellm/work/src/hpstr/build/processors && /opt/rh/devtoolset-8/root/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /u/re/alspellm/work/src/hpstr/processors/src/TrackingProcessor.cxx > CMakeFiles/processors.dir/src/TrackingProcessor.cxx.i

processors/CMakeFiles/processors.dir/src/TrackingProcessor.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/processors.dir/src/TrackingProcessor.cxx.s"
	cd /u/re/alspellm/work/src/hpstr/build/processors && /opt/rh/devtoolset-8/root/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /u/re/alspellm/work/src/hpstr/processors/src/TrackingProcessor.cxx -o CMakeFiles/processors.dir/src/TrackingProcessor.cxx.s

processors/CMakeFiles/processors.dir/src/utilities.cxx.o: processors/CMakeFiles/processors.dir/flags.make
processors/CMakeFiles/processors.dir/src/utilities.cxx.o: ../processors/src/utilities.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/u/re/alspellm/work/src/hpstr/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building CXX object processors/CMakeFiles/processors.dir/src/utilities.cxx.o"
	cd /u/re/alspellm/work/src/hpstr/build/processors && /opt/rh/devtoolset-8/root/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/processors.dir/src/utilities.cxx.o -c /u/re/alspellm/work/src/hpstr/processors/src/utilities.cxx

processors/CMakeFiles/processors.dir/src/utilities.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/processors.dir/src/utilities.cxx.i"
	cd /u/re/alspellm/work/src/hpstr/build/processors && /opt/rh/devtoolset-8/root/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /u/re/alspellm/work/src/hpstr/processors/src/utilities.cxx > CMakeFiles/processors.dir/src/utilities.cxx.i

processors/CMakeFiles/processors.dir/src/utilities.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/processors.dir/src/utilities.cxx.s"
	cd /u/re/alspellm/work/src/hpstr/build/processors && /opt/rh/devtoolset-8/root/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /u/re/alspellm/work/src/hpstr/processors/src/utilities.cxx -o CMakeFiles/processors.dir/src/utilities.cxx.s

# Object files for target processors
processors_OBJECTS = \
"CMakeFiles/processors.dir/src/ClusterOnTrackAnaProcessor.cxx.o" \
"CMakeFiles/processors.dir/src/ECalDataProcessor.cxx.o" \
"CMakeFiles/processors.dir/src/EventProcessor.cxx.o" \
"CMakeFiles/processors.dir/src/HPSEventProcessor.cxx.o" \
"CMakeFiles/processors.dir/src/ParticleProcessor.cxx.o" \
"CMakeFiles/processors.dir/src/RefittedTracksProcessor.cxx.o" \
"CMakeFiles/processors.dir/src/SvtCondAnaProcessor.cxx.o" \
"CMakeFiles/processors.dir/src/SvtDataProcessor.cxx.o" \
"CMakeFiles/processors.dir/src/SvtRawDataProcessor.cxx.o" \
"CMakeFiles/processors.dir/src/TrackingProcessor.cxx.o" \
"CMakeFiles/processors.dir/src/utilities.cxx.o"

# External object files for target processors
processors_EXTERNAL_OBJECTS =

processors/libprocessors.so: processors/CMakeFiles/processors.dir/src/ClusterOnTrackAnaProcessor.cxx.o
processors/libprocessors.so: processors/CMakeFiles/processors.dir/src/ECalDataProcessor.cxx.o
processors/libprocessors.so: processors/CMakeFiles/processors.dir/src/EventProcessor.cxx.o
processors/libprocessors.so: processors/CMakeFiles/processors.dir/src/HPSEventProcessor.cxx.o
processors/libprocessors.so: processors/CMakeFiles/processors.dir/src/ParticleProcessor.cxx.o
processors/libprocessors.so: processors/CMakeFiles/processors.dir/src/RefittedTracksProcessor.cxx.o
processors/libprocessors.so: processors/CMakeFiles/processors.dir/src/SvtCondAnaProcessor.cxx.o
processors/libprocessors.so: processors/CMakeFiles/processors.dir/src/SvtDataProcessor.cxx.o
processors/libprocessors.so: processors/CMakeFiles/processors.dir/src/SvtRawDataProcessor.cxx.o
processors/libprocessors.so: processors/CMakeFiles/processors.dir/src/TrackingProcessor.cxx.o
processors/libprocessors.so: processors/CMakeFiles/processors.dir/src/utilities.cxx.o
processors/libprocessors.so: processors/CMakeFiles/processors.dir/build.make
processors/libprocessors.so: processing/libprocessing.so
processors/libprocessors.so: analysis/libanalysis.so
processors/libprocessors.so: /nfs/slac/g/hps3/users/bravo/src/root/buildV61204/lib/libCore.so
processors/libprocessors.so: /nfs/slac/g/hps3/users/bravo/src/root/buildV61204/lib/libImt.so
processors/libprocessors.so: /nfs/slac/g/hps3/users/bravo/src/root/buildV61204/lib/libRIO.so
processors/libprocessors.so: /nfs/slac/g/hps3/users/bravo/src/root/buildV61204/lib/libNet.so
processors/libprocessors.so: /nfs/slac/g/hps3/users/bravo/src/root/buildV61204/lib/libHist.so
processors/libprocessors.so: /nfs/slac/g/hps3/users/bravo/src/root/buildV61204/lib/libGraf.so
processors/libprocessors.so: /nfs/slac/g/hps3/users/bravo/src/root/buildV61204/lib/libGraf3d.so
processors/libprocessors.so: /nfs/slac/g/hps3/users/bravo/src/root/buildV61204/lib/libGpad.so
processors/libprocessors.so: /nfs/slac/g/hps3/users/bravo/src/root/buildV61204/lib/libTree.so
processors/libprocessors.so: /nfs/slac/g/hps3/users/bravo/src/root/buildV61204/lib/libTreePlayer.so
processors/libprocessors.so: /nfs/slac/g/hps3/users/bravo/src/root/buildV61204/lib/libRint.so
processors/libprocessors.so: /nfs/slac/g/hps3/users/bravo/src/root/buildV61204/lib/libPostscript.so
processors/libprocessors.so: /nfs/slac/g/hps3/users/bravo/src/root/buildV61204/lib/libMatrix.so
processors/libprocessors.so: /nfs/slac/g/hps3/users/bravo/src/root/buildV61204/lib/libPhysics.so
processors/libprocessors.so: /nfs/slac/g/hps3/users/bravo/src/root/buildV61204/lib/libMathCore.so
processors/libprocessors.so: /nfs/slac/g/hps3/users/bravo/src/root/buildV61204/lib/libThread.so
processors/libprocessors.so: /nfs/slac/g/hps3/users/bravo/src/root/buildV61204/lib/libMultiProc.so
processors/libprocessors.so: /nfs/slac/g/hps3/users/bravo/src/root/buildV61204/lib/libPyROOT.so
processors/libprocessors.so: /nfs/slac/g/hps3/users/bravo/src/root/buildV61204/lib/libGeom.so
processors/libprocessors.so: /nfs/slac/g/hps3/users/bravo/src/root/buildV61204/lib/libEve.so
processors/libprocessors.so: /nfs/slac/g/hps3/users/bravo/src/root/buildV61204/lib/libGui.so
processors/libprocessors.so: /nfs/slac/g/hps3/users/alspellm/src/LCIO/install/lib/liblcio.so
processors/libprocessors.so: /usr/lib64/libpython2.7.so
processors/libprocessors.so: event/libevent.so
processors/libprocessors.so: /nfs/slac/g/hps3/users/alspellm/src/LCIO/install/lib/liblcio.so
processors/libprocessors.so: /nfs/slac/g/hps3/users/bravo/src/root/buildV61204/lib/libCore.so
processors/libprocessors.so: /nfs/slac/g/hps3/users/bravo/src/root/buildV61204/lib/libImt.so
processors/libprocessors.so: /nfs/slac/g/hps3/users/bravo/src/root/buildV61204/lib/libRIO.so
processors/libprocessors.so: /nfs/slac/g/hps3/users/bravo/src/root/buildV61204/lib/libNet.so
processors/libprocessors.so: /nfs/slac/g/hps3/users/bravo/src/root/buildV61204/lib/libHist.so
processors/libprocessors.so: /nfs/slac/g/hps3/users/bravo/src/root/buildV61204/lib/libGraf.so
processors/libprocessors.so: /nfs/slac/g/hps3/users/bravo/src/root/buildV61204/lib/libGraf3d.so
processors/libprocessors.so: /nfs/slac/g/hps3/users/bravo/src/root/buildV61204/lib/libGpad.so
processors/libprocessors.so: /nfs/slac/g/hps3/users/bravo/src/root/buildV61204/lib/libTree.so
processors/libprocessors.so: /nfs/slac/g/hps3/users/bravo/src/root/buildV61204/lib/libTreePlayer.so
processors/libprocessors.so: /nfs/slac/g/hps3/users/bravo/src/root/buildV61204/lib/libRint.so
processors/libprocessors.so: /nfs/slac/g/hps3/users/bravo/src/root/buildV61204/lib/libPostscript.so
processors/libprocessors.so: /nfs/slac/g/hps3/users/bravo/src/root/buildV61204/lib/libMatrix.so
processors/libprocessors.so: /nfs/slac/g/hps3/users/bravo/src/root/buildV61204/lib/libPhysics.so
processors/libprocessors.so: /nfs/slac/g/hps3/users/bravo/src/root/buildV61204/lib/libMathCore.so
processors/libprocessors.so: /nfs/slac/g/hps3/users/bravo/src/root/buildV61204/lib/libThread.so
processors/libprocessors.so: /nfs/slac/g/hps3/users/bravo/src/root/buildV61204/lib/libMultiProc.so
processors/libprocessors.so: /nfs/slac/g/hps3/users/bravo/src/root/buildV61204/lib/libPyROOT.so
processors/libprocessors.so: /nfs/slac/g/hps3/users/bravo/src/root/buildV61204/lib/libGeom.so
processors/libprocessors.so: /nfs/slac/g/hps3/users/bravo/src/root/buildV61204/lib/libEve.so
processors/libprocessors.so: /nfs/slac/g/hps3/users/bravo/src/root/buildV61204/lib/libGui.so
processors/libprocessors.so: processors/CMakeFiles/processors.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/u/re/alspellm/work/src/hpstr/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Linking CXX shared library libprocessors.so"
	cd /u/re/alspellm/work/src/hpstr/build/processors && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/processors.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
processors/CMakeFiles/processors.dir/build: processors/libprocessors.so

.PHONY : processors/CMakeFiles/processors.dir/build

processors/CMakeFiles/processors.dir/clean:
	cd /u/re/alspellm/work/src/hpstr/build/processors && $(CMAKE_COMMAND) -P CMakeFiles/processors.dir/cmake_clean.cmake
.PHONY : processors/CMakeFiles/processors.dir/clean

processors/CMakeFiles/processors.dir/depend:
	cd /u/re/alspellm/work/src/hpstr/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /u/re/alspellm/work/src/hpstr /u/re/alspellm/work/src/hpstr/processors /u/re/alspellm/work/src/hpstr/build /u/re/alspellm/work/src/hpstr/build/processors /u/re/alspellm/work/src/hpstr/build/processors/CMakeFiles/processors.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : processors/CMakeFiles/processors.dir/depend

