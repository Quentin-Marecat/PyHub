cmake_minimum_required (VERSION 3.5)
project ("PyHub")
enable_language(Fortran)

# Currently setting the Fortran compiler to use -std=gnu, change this if you
# want a specific standard
set(FVERSION "-std=gnu")

# BLAS
if (NOT BLAS_LIBRARIES)
find_package(BLAS)
endif()
if (NOT BLAS_LIBRARIES)
    message(FATAL_ERROR "A required library with BLAS API not found.")
else()
    message(STATUS "BLAS libraries: ${BLAS_LIBRARIES}")
endif()

if (APPLE)
    include_directories("${BLAS_LIBRARIES}/Frameworks/vecLib.framework/Headers")
endif()

if (APPLE)
    set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
    set(CMAKE_INSTALL_RPATH "@loader_path;@loader_path/deps/lib;@loader_path/deps/lib64")
    set(CMAKE_BUILD_RPATH "@loader_path;@loader_path/deps/lib;@loader_path/deps/lib64")
else()
    set(CMAKE_SKIP_BUILD_RPATH  True)
    set(CMAKE_BUILD_WITH_INSTALL_RPATH True)
    set(CMAKE_INSTALL_RPATH "\$ORIGIN:\$ORIGIN/deps/lib:\$ORIGIN/deps/lib64")
endif()

# HDF5
find_package(HDF5 REQUIRED COMPONENTS Fortran HL)
if (NOT HDF5_LIBRARIES)
    message(FATAL_ERROR "A required library with HDF5 API not found.")
else()
    message(STATUS "HDF5 libraries: ${HDF5_LIBRARIES}")
endif()


# Fortran FLAGS
set(CMAKE_Fortran_FLAGS "-O2 -g3 -Wall -Wextra -Wconversion -fbounds-check -fbacktrace -fcheck=all ${CMAKE_Fortran_FLAGS} ${FVERSION}")
include(CheckFortranCompilerFlag)
CHECK_Fortran_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
if (COMPILER_SUPPORTS_MARCH_NATIVE)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -march=native")
endif()

add_subdirectory(src)

#install(TARGETS hubbard.x RUNTIME DESTINATION "bin")
