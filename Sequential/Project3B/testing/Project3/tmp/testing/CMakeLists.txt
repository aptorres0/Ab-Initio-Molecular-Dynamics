cmake_minimum_required(VERSION 3.10.1)

#=== Set name of our project
project(project3)

set(CMAKE_CXX_STANDARD 17)

#=== Find all source files in a directory
aux_source_directory(./src SOURCE_FILES)

set(CMAKE_CXX_FLAGS "-Wall -Werror -pedantic -Wextra -O3")

# Bring headers into the project
include_directories(include)

# Add Sources
file( GLOB SOURCES "src/*")

set( LIBRARY_OUTPUT_PATH "${CMAKE_CURRENT_SOURCE_DIR}/build/lib" )
add_library( project3 SHARED ${SOURCES})

set( EXECUTABLE_OUTPUT_PATH "${CMAKE_CURRENT_SOURCE_DIR}/build/scripts")
add_executable( main scripts/main.cxx )
target_link_libraries(main project3)

