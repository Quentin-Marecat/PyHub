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
CMAKE_SOURCE_DIR = /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs

# Include any dependencies generated for this target.
include src/CMakeFiles/hubbard.x.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include src/CMakeFiles/hubbard.x.dir/compiler_depend.make

# Include the progress variables for this target.
include src/CMakeFiles/hubbard.x.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/hubbard.x.dir/flags.make

src/CMakeFiles/hubbard.x.dir/hubmod.f90.o: src/CMakeFiles/hubbard.x.dir/flags.make
src/CMakeFiles/hubbard.x.dir/hubmod.f90.o: src/hubmod.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object src/CMakeFiles/hubbard.x.dir/hubmod.f90.o"
	cd /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src/hubmod.f90 -o CMakeFiles/hubbard.x.dir/hubmod.f90.o

src/CMakeFiles/hubbard.x.dir/hubmod.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/hubbard.x.dir/hubmod.f90.i"
	cd /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src/hubmod.f90 > CMakeFiles/hubbard.x.dir/hubmod.f90.i

src/CMakeFiles/hubbard.x.dir/hubmod.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/hubbard.x.dir/hubmod.f90.s"
	cd /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src/hubmod.f90 -o CMakeFiles/hubbard.x.dir/hubmod.f90.s

src/CMakeFiles/hubbard.x.dir/funcmod.f90.o: src/CMakeFiles/hubbard.x.dir/flags.make
src/CMakeFiles/hubbard.x.dir/funcmod.f90.o: src/funcmod.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building Fortran object src/CMakeFiles/hubbard.x.dir/funcmod.f90.o"
	cd /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src/funcmod.f90 -o CMakeFiles/hubbard.x.dir/funcmod.f90.o

src/CMakeFiles/hubbard.x.dir/funcmod.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/hubbard.x.dir/funcmod.f90.i"
	cd /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src/funcmod.f90 > CMakeFiles/hubbard.x.dir/funcmod.f90.i

src/CMakeFiles/hubbard.x.dir/funcmod.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/hubbard.x.dir/funcmod.f90.s"
	cd /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src/funcmod.f90 -o CMakeFiles/hubbard.x.dir/funcmod.f90.s

src/CMakeFiles/hubbard.x.dir/basismod.f90.o: src/CMakeFiles/hubbard.x.dir/flags.make
src/CMakeFiles/hubbard.x.dir/basismod.f90.o: src/basismod.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building Fortran object src/CMakeFiles/hubbard.x.dir/basismod.f90.o"
	cd /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src/basismod.f90 -o CMakeFiles/hubbard.x.dir/basismod.f90.o

src/CMakeFiles/hubbard.x.dir/basismod.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/hubbard.x.dir/basismod.f90.i"
	cd /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src/basismod.f90 > CMakeFiles/hubbard.x.dir/basismod.f90.i

src/CMakeFiles/hubbard.x.dir/basismod.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/hubbard.x.dir/basismod.f90.s"
	cd /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src/basismod.f90 -o CMakeFiles/hubbard.x.dir/basismod.f90.s

src/CMakeFiles/hubbard.x.dir/hubbard.f90.o: src/CMakeFiles/hubbard.x.dir/flags.make
src/CMakeFiles/hubbard.x.dir/hubbard.f90.o: src/hubbard.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building Fortran object src/CMakeFiles/hubbard.x.dir/hubbard.f90.o"
	cd /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src/hubbard.f90 -o CMakeFiles/hubbard.x.dir/hubbard.f90.o

src/CMakeFiles/hubbard.x.dir/hubbard.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/hubbard.x.dir/hubbard.f90.i"
	cd /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src/hubbard.f90 > CMakeFiles/hubbard.x.dir/hubbard.f90.i

src/CMakeFiles/hubbard.x.dir/hubbard.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/hubbard.x.dir/hubbard.f90.s"
	cd /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src/hubbard.f90 -o CMakeFiles/hubbard.x.dir/hubbard.f90.s

src/CMakeFiles/hubbard.x.dir/diagmat.f90.o: src/CMakeFiles/hubbard.x.dir/flags.make
src/CMakeFiles/hubbard.x.dir/diagmat.f90.o: src/diagmat.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building Fortran object src/CMakeFiles/hubbard.x.dir/diagmat.f90.o"
	cd /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src/diagmat.f90 -o CMakeFiles/hubbard.x.dir/diagmat.f90.o

src/CMakeFiles/hubbard.x.dir/diagmat.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/hubbard.x.dir/diagmat.f90.i"
	cd /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src/diagmat.f90 > CMakeFiles/hubbard.x.dir/diagmat.f90.i

src/CMakeFiles/hubbard.x.dir/diagmat.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/hubbard.x.dir/diagmat.f90.s"
	cd /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src/diagmat.f90 -o CMakeFiles/hubbard.x.dir/diagmat.f90.s

src/CMakeFiles/hubbard.x.dir/basis.f90.o: src/CMakeFiles/hubbard.x.dir/flags.make
src/CMakeFiles/hubbard.x.dir/basis.f90.o: src/basis.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building Fortran object src/CMakeFiles/hubbard.x.dir/basis.f90.o"
	cd /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src/basis.f90 -o CMakeFiles/hubbard.x.dir/basis.f90.o

src/CMakeFiles/hubbard.x.dir/basis.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/hubbard.x.dir/basis.f90.i"
	cd /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src/basis.f90 > CMakeFiles/hubbard.x.dir/basis.f90.i

src/CMakeFiles/hubbard.x.dir/basis.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/hubbard.x.dir/basis.f90.s"
	cd /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src/basis.f90 -o CMakeFiles/hubbard.x.dir/basis.f90.s

src/CMakeFiles/hubbard.x.dir/hprod.f90.o: src/CMakeFiles/hubbard.x.dir/flags.make
src/CMakeFiles/hubbard.x.dir/hprod.f90.o: src/hprod.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building Fortran object src/CMakeFiles/hubbard.x.dir/hprod.f90.o"
	cd /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src/hprod.f90 -o CMakeFiles/hubbard.x.dir/hprod.f90.o

src/CMakeFiles/hubbard.x.dir/hprod.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/hubbard.x.dir/hprod.f90.i"
	cd /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src/hprod.f90 > CMakeFiles/hubbard.x.dir/hprod.f90.i

src/CMakeFiles/hubbard.x.dir/hprod.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/hubbard.x.dir/hprod.f90.s"
	cd /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src/hprod.f90 -o CMakeFiles/hubbard.x.dir/hprod.f90.s

src/CMakeFiles/hubbard.x.dir/lanczos.f90.o: src/CMakeFiles/hubbard.x.dir/flags.make
src/CMakeFiles/hubbard.x.dir/lanczos.f90.o: src/lanczos.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building Fortran object src/CMakeFiles/hubbard.x.dir/lanczos.f90.o"
	cd /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src/lanczos.f90 -o CMakeFiles/hubbard.x.dir/lanczos.f90.o

src/CMakeFiles/hubbard.x.dir/lanczos.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/hubbard.x.dir/lanczos.f90.i"
	cd /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src/lanczos.f90 > CMakeFiles/hubbard.x.dir/lanczos.f90.i

src/CMakeFiles/hubbard.x.dir/lanczos.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/hubbard.x.dir/lanczos.f90.s"
	cd /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src/lanczos.f90 -o CMakeFiles/hubbard.x.dir/lanczos.f90.s

src/CMakeFiles/hubbard.x.dir/solve_hubbard.f90.o: src/CMakeFiles/hubbard.x.dir/flags.make
src/CMakeFiles/hubbard.x.dir/solve_hubbard.f90.o: src/solve_hubbard.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building Fortran object src/CMakeFiles/hubbard.x.dir/solve_hubbard.f90.o"
	cd /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src/solve_hubbard.f90 -o CMakeFiles/hubbard.x.dir/solve_hubbard.f90.o

src/CMakeFiles/hubbard.x.dir/solve_hubbard.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/hubbard.x.dir/solve_hubbard.f90.i"
	cd /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src/solve_hubbard.f90 > CMakeFiles/hubbard.x.dir/solve_hubbard.f90.i

src/CMakeFiles/hubbard.x.dir/solve_hubbard.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/hubbard.x.dir/solve_hubbard.f90.s"
	cd /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src/solve_hubbard.f90 -o CMakeFiles/hubbard.x.dir/solve_hubbard.f90.s

src/CMakeFiles/hubbard.x.dir/static_rq.f90.o: src/CMakeFiles/hubbard.x.dir/flags.make
src/CMakeFiles/hubbard.x.dir/static_rq.f90.o: src/static_rq.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building Fortran object src/CMakeFiles/hubbard.x.dir/static_rq.f90.o"
	cd /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src/static_rq.f90 -o CMakeFiles/hubbard.x.dir/static_rq.f90.o

src/CMakeFiles/hubbard.x.dir/static_rq.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/hubbard.x.dir/static_rq.f90.i"
	cd /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src/static_rq.f90 > CMakeFiles/hubbard.x.dir/static_rq.f90.i

src/CMakeFiles/hubbard.x.dir/static_rq.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/hubbard.x.dir/static_rq.f90.s"
	cd /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src/static_rq.f90 -o CMakeFiles/hubbard.x.dir/static_rq.f90.s

src/CMakeFiles/hubbard.x.dir/exact_diag.f90.o: src/CMakeFiles/hubbard.x.dir/flags.make
src/CMakeFiles/hubbard.x.dir/exact_diag.f90.o: src/exact_diag.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building Fortran object src/CMakeFiles/hubbard.x.dir/exact_diag.f90.o"
	cd /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src/exact_diag.f90 -o CMakeFiles/hubbard.x.dir/exact_diag.f90.o

src/CMakeFiles/hubbard.x.dir/exact_diag.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/hubbard.x.dir/exact_diag.f90.i"
	cd /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src/exact_diag.f90 > CMakeFiles/hubbard.x.dir/exact_diag.f90.i

src/CMakeFiles/hubbard.x.dir/exact_diag.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/hubbard.x.dir/exact_diag.f90.s"
	cd /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src/exact_diag.f90 -o CMakeFiles/hubbard.x.dir/exact_diag.f90.s

src/CMakeFiles/hubbard.x.dir/boltzmann.f90.o: src/CMakeFiles/hubbard.x.dir/flags.make
src/CMakeFiles/hubbard.x.dir/boltzmann.f90.o: src/boltzmann.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Building Fortran object src/CMakeFiles/hubbard.x.dir/boltzmann.f90.o"
	cd /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src/boltzmann.f90 -o CMakeFiles/hubbard.x.dir/boltzmann.f90.o

src/CMakeFiles/hubbard.x.dir/boltzmann.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/hubbard.x.dir/boltzmann.f90.i"
	cd /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src/boltzmann.f90 > CMakeFiles/hubbard.x.dir/boltzmann.f90.i

src/CMakeFiles/hubbard.x.dir/boltzmann.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/hubbard.x.dir/boltzmann.f90.s"
	cd /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src/boltzmann.f90 -o CMakeFiles/hubbard.x.dir/boltzmann.f90.s

src/CMakeFiles/hubbard.x.dir/spgf_ed.f90.o: src/CMakeFiles/hubbard.x.dir/flags.make
src/CMakeFiles/hubbard.x.dir/spgf_ed.f90.o: src/spgf_ed.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/CMakeFiles --progress-num=$(CMAKE_PROGRESS_13) "Building Fortran object src/CMakeFiles/hubbard.x.dir/spgf_ed.f90.o"
	cd /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src/spgf_ed.f90 -o CMakeFiles/hubbard.x.dir/spgf_ed.f90.o

src/CMakeFiles/hubbard.x.dir/spgf_ed.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/hubbard.x.dir/spgf_ed.f90.i"
	cd /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src/spgf_ed.f90 > CMakeFiles/hubbard.x.dir/spgf_ed.f90.i

src/CMakeFiles/hubbard.x.dir/spgf_ed.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/hubbard.x.dir/spgf_ed.f90.s"
	cd /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src/spgf_ed.f90 -o CMakeFiles/hubbard.x.dir/spgf_ed.f90.s

src/CMakeFiles/hubbard.x.dir/spgf_bandlanczos.f90.o: src/CMakeFiles/hubbard.x.dir/flags.make
src/CMakeFiles/hubbard.x.dir/spgf_bandlanczos.f90.o: src/spgf_bandlanczos.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/CMakeFiles --progress-num=$(CMAKE_PROGRESS_14) "Building Fortran object src/CMakeFiles/hubbard.x.dir/spgf_bandlanczos.f90.o"
	cd /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src/spgf_bandlanczos.f90 -o CMakeFiles/hubbard.x.dir/spgf_bandlanczos.f90.o

src/CMakeFiles/hubbard.x.dir/spgf_bandlanczos.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/hubbard.x.dir/spgf_bandlanczos.f90.i"
	cd /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src/spgf_bandlanczos.f90 > CMakeFiles/hubbard.x.dir/spgf_bandlanczos.f90.i

src/CMakeFiles/hubbard.x.dir/spgf_bandlanczos.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/hubbard.x.dir/spgf_bandlanczos.f90.s"
	cd /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src/spgf_bandlanczos.f90 -o CMakeFiles/hubbard.x.dir/spgf_bandlanczos.f90.s

# Object files for target hubbard.x
hubbard_x_OBJECTS = \
"CMakeFiles/hubbard.x.dir/hubmod.f90.o" \
"CMakeFiles/hubbard.x.dir/funcmod.f90.o" \
"CMakeFiles/hubbard.x.dir/basismod.f90.o" \
"CMakeFiles/hubbard.x.dir/hubbard.f90.o" \
"CMakeFiles/hubbard.x.dir/diagmat.f90.o" \
"CMakeFiles/hubbard.x.dir/basis.f90.o" \
"CMakeFiles/hubbard.x.dir/hprod.f90.o" \
"CMakeFiles/hubbard.x.dir/lanczos.f90.o" \
"CMakeFiles/hubbard.x.dir/solve_hubbard.f90.o" \
"CMakeFiles/hubbard.x.dir/static_rq.f90.o" \
"CMakeFiles/hubbard.x.dir/exact_diag.f90.o" \
"CMakeFiles/hubbard.x.dir/boltzmann.f90.o" \
"CMakeFiles/hubbard.x.dir/spgf_ed.f90.o" \
"CMakeFiles/hubbard.x.dir/spgf_bandlanczos.f90.o"

# External object files for target hubbard.x
hubbard_x_EXTERNAL_OBJECTS =

src/hubbard.x: src/CMakeFiles/hubbard.x.dir/hubmod.f90.o
src/hubbard.x: src/CMakeFiles/hubbard.x.dir/funcmod.f90.o
src/hubbard.x: src/CMakeFiles/hubbard.x.dir/basismod.f90.o
src/hubbard.x: src/CMakeFiles/hubbard.x.dir/hubbard.f90.o
src/hubbard.x: src/CMakeFiles/hubbard.x.dir/diagmat.f90.o
src/hubbard.x: src/CMakeFiles/hubbard.x.dir/basis.f90.o
src/hubbard.x: src/CMakeFiles/hubbard.x.dir/hprod.f90.o
src/hubbard.x: src/CMakeFiles/hubbard.x.dir/lanczos.f90.o
src/hubbard.x: src/CMakeFiles/hubbard.x.dir/solve_hubbard.f90.o
src/hubbard.x: src/CMakeFiles/hubbard.x.dir/static_rq.f90.o
src/hubbard.x: src/CMakeFiles/hubbard.x.dir/exact_diag.f90.o
src/hubbard.x: src/CMakeFiles/hubbard.x.dir/boltzmann.f90.o
src/hubbard.x: src/CMakeFiles/hubbard.x.dir/spgf_ed.f90.o
src/hubbard.x: src/CMakeFiles/hubbard.x.dir/spgf_bandlanczos.f90.o
src/hubbard.x: src/CMakeFiles/hubbard.x.dir/build.make
src/hubbard.x: /usr/lib/x86_64-linux-gnu/libopenblas.so
src/hubbard.x: /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_fortran.so
src/hubbard.x: src/CMakeFiles/hubbard.x.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/CMakeFiles --progress-num=$(CMAKE_PROGRESS_15) "Linking Fortran executable hubbard.x"
	cd /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/hubbard.x.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/hubbard.x.dir/build: src/hubbard.x
.PHONY : src/CMakeFiles/hubbard.x.dir/build

src/CMakeFiles/hubbard.x.dir/clean:
	cd /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src && $(CMAKE_COMMAND) -P CMakeFiles/hubbard.x.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/hubbard.x.dir/clean

src/CMakeFiles/hubbard.x.dir/depend:
	cd /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src /home/qmarecat/Desktop/Work/Github/PyHub/pyhub/libs/src/CMakeFiles/hubbard.x.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/hubbard.x.dir/depend

