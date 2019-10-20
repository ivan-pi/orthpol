cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

enable_language(Fortran)

project(oc)

#set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
# Set output paths for modules, archives, and executables
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Make sure that CMAKE_BUILD_TYPE is either Debug or Release:
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug
        CACHE STRING "Build type (Debug, Release)" FORCE)
endif ()
if (NOT (CMAKE_BUILD_TYPE STREQUAL "Debug" OR
        CMAKE_BUILD_TYPE STREQUAL "Release"))
    message("${CMAKE_BUILD_TYPE}")
    message(FATAL_ERROR "CMAKE_BUILD_TYPE must be one of: Debug, Release (current value: '${CMAKE_BUILD_TYPE}')")
endif ()

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    execute_process(COMMAND ${CMAKE_Fortran_COMPILER} -dumpversion
                        OUTPUT_VARIABLE GCC_VERSION)
    #if (GCC_VERSION VERSION_LESS 4.5)
        set(common "-std=f2008 -Wall -Wextra -Wimplicit-interface -fPIC -fmax-errors=1")
        set(CMAKE_Fortran_FLAGS_RELEASE "${common} -O3 -march=native -ffast-math -funroll-loops")
        set(CMAKE_Fortran_FLAGS_DEBUG   "${common} -g -fbounds-check -fcheck-array-temporaries -fbacktrace")
    #endif()

    # gfortran
    # Enable this if you want to check for single/double corruption (and use
    # the Debug build):
    #set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -fdefault-real-8")
endif ()

enable_testing()

add_subdirectory(src)
add_subdirectory(tests)



message("\n")
message("Configuration results")
message("---------------------")
message("Fortran compiler: ${CMAKE_Fortran_COMPILER}")
message("Build type: ${CMAKE_BUILD_TYPE}")
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    message("Fortran compiler flags: ${CMAKE_Fortran_FLAGS_DEBUG}")
else ()
    message("Fortran compiler flags: ${CMAKE_Fortran_FLAGS_RELEASE}")
endif ()
message("Installation prefix: ${CMAKE_INSTALL_PREFIX}")