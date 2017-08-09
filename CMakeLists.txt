# Original file Copyright CMU462 Fall 2015: 
# Kayvon Fatahalian, Keenan Crane,
# Sky Gao, Bryce Summers, Michael Choquette.
cmake_minimum_required(VERSION 2.8)
project(sph)

#-------------------------------------------------------------------------------
# Build options
#-------------------------------------------------------------------------------
option(BUILD_DEBUG     "Build with debug settings"    OFF)
option(BUILD_DOCS      "Build documentation"          OFF)

#-------------------------------------------------------------------------------
# Platform-specific settings
#-------------------------------------------------------------------------------

###################
# Building on OSX #
###################
if(APPLE)

  # Clang configuration
  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")

    set(CLANG_CXX_FLAGS "-std=c++11 -m64")

    if(BUILD_DEBUG)
      set(CMAKE_BUILD_TYPE Debug)
    else(BUILD_DEBUG)
      set(CLANG_CXX_FLAGS "${CLANG_CXX_FLAGS} -O3")
      set(CLANG_CXX_FLAGS "${CLANG_CXX_FLAGS} -funroll-loops")
      set(CLANG_CXX_FLAGS "${CLANG_CXX_FLAGS} -Wno-narrowing")
      set(CLANG_CXX_FLAGS "${CLANG_CXX_FLAGS} -Wno-deprecated-register")
    endif(BUILD_DEBUG)

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CLANG_CXX_FLAGS}")

  endif()

  # GCC configuration
  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")

    set(GCC_CXX_FLAGS "-std=c++11 -m64")

    if(BUILD_DEBUG)
      set(CMAKE_BUILD_TYPE Debug)
    else(BUILD_DEBUG)
      set(GCC_CXX_FLAGS "${GCC_CXX_FLAGS} -O3")
      set(GCC_CXX_FLAGS "${GCC_CXX_FLAGS} -fopenmp")
    endif(BUILD_DEBUG)

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_CXX_FLAGS}")

  endif()

endif(APPLE)

##################
# Build on Linux #
##################
set(LINUX UNIX AND NOT APPLE)
if(LINUX)

  # GCC only
  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")

    set(GCC_CXX_FLAGS "-std=c++11 -m64")

    # Debug configuration
    if(BUILD_DEBUG)
        set(CMAKE_BUILD_TYPE Debug)
    else(BUILD_DEBUG)
        set(GCC_CXX_FLAGS "${GCC_CXX_FLAGS} -O3")
        set(GCC_CXX_FLAGS "${GCC_CXX_FLAGS} -fopenmp")
    endif(BUILD_DEBUG)

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_CXX_FLAGS}")

  endif()

endif(LINUX)

####################
# Build on Windows #
####################
if(WIN32)

  if(MSVC)

    set(MSVC_CXX_FLAGS "-std=c++11")

    if(BUILD_DEBUG)
        set(CMAKE_BUILD_TYPE Debug)
    else(BUILD_DEBUG)
    endif(BUILD_DEBUG)

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MSVC_CXX_FLAGS}")

  endif(MSVC)

  if(MINGW)

    set(MSVC_CXX_FLAGS "-std=c++11")

    if(BUILD_DEBUG)
        set(CMAKE_BUILD_TYPE Debug)
    else(BUILD_DEBUG)
    endif(BUILD_DEBUG)

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MSVC_CXX_FLAGS}")

  endif(MINGW)

endif(WIN32)

#-------------------------------------------------------------------------------
# Find dependencies
#-------------------------------------------------------------------------------

# Required packages
find_package(OpenGL REQUIRED)
find_package(Threads REQUIRED)
find_package(Freetype REQUIRED)

# glfw
add_subdirectory(glfw-3.2.1)
include_directories(glfw-3.2.1/include)

# glew
add_subdirectory(glew)
include_directories(glew/include)
# set_property( TARGET glew APPEND_STRING PROPERTY COMPILE_FLAGS -w )

#-------------------------------------------------------------------------------
# Add subdirectories
#-------------------------------------------------------------------------------
add_subdirectory(src)

# build documentation
#if(BUILD_DOCS)
#  find_package(DOXYGEN)
#  if(DOXYGEN_FOUND AND BUILD_DOCS)
#    add_subdirectory(docs)
#  endif()
#endif()

# Install settings
set(CMAKE_INSTALL_PREFIX “${sph_SOURCE_DIR}/“)
