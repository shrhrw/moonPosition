# Minimum version of CMake required to build this project
cmake_minimum_required( VERSION 3.22.1 )

# Name of project
project( Moon LANGUAGES CXX C )

# CMake standard (C++11)
set( CMAKE_CXX_STANDARD 11 )

# Retrieve googletest package
include( FetchContent )
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG main
)

FetchContent_MakeAvailable( googletest )

# Enable googletest
enable_testing()

# Add path to all dependency directories
include_directories( 
    include
    src
    test
)

# Add all the source files needed to build the executable
add_library( JulianDate STATIC src/JD.c )
add_library( MoonPositionModel STATIC src/moonPosition.cpp )

add_executable( MoonPositionTests test/src/moonPositionTests.cpp )

target_link_libraries( MoonPositionTests
                       MoonPositionModel
                       JulianDate
                       gtest_main
                     )

include( GoogleTest )
gtest_discover_tests( MoonPositionTests )
