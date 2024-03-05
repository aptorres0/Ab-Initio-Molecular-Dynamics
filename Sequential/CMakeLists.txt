cmake_minimum_required(VERSION 3.10.1)

#=== Set name of our project
project(project3)

#=== Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT 6)

#=== Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_STANDARD 20)

#find_package(OpenMP REQUIRED)
#set(CMAKE_CXX_FLAGS "-Wall -Werror -Wextra -Wpedantic -g -Og")
#set(CMAKE_CXX_FLAGS "-Wall -Werror -Wextra -Wpedantic -O0 -g")
#set(CMAKE_CXX_FLAGS "-Wall -Werror -pedantic -Wextra -Ofast -funroll-loops")
set(CMAKE_CXX_FLAGS "-Wall -Werror -pedantic -Wextra -O3")
#set(CMAKE_CXX_FLAGS "-O3")

# See: https://stackoverflow.com/questions/14492436/g-optimization-beyond-o3-ofast
#set(CMAKE_CXX_FLAGS "-Wall -Werror -pedantic -Wextra -Ofast -mcpu=native -flto -fno-signed-zeros -fno-trapping-math -funroll-loops")

# Bring headers into the project
include_directories(include)

# Add Sources
file( GLOB SOURCES "src/*")

set( LIBRARY_OUTPUT_PATH "${CMAKE_CURRENT_SOURCE_DIR}/build/lib" )
add_library( project3 SHARED ${SOURCES})
#
target_link_libraries(project3 ${ROOT_LIBRARIES})

set( EXECUTABLE_OUTPUT_PATH "${CMAKE_CURRENT_SOURCE_DIR}/build/scripts")
add_executable( main scripts/main.cxx )
target_link_libraries(main project3 ${ROOT_LIBRARIES})
#target_link_libraries(main project3 OpenMP::OpenMP_CXX)

#add_executable( newmain scripts/newmain.cxx )
#target_link_libraries(newmain project3 ${ROOT_LIBRARIES})
#target_link_libraries(newmain project3 OpenMP::OpenMP_CXX)

#add_executable( test scripts/test.cxx )
#target_link_libraries(test project3)
#target_link_libraries(test project3 ${ROOT_LIBRARIES})

message("CMAKE_BUILD_TYPE is ${CMAKE_BUILD_TYPE}")
message("CMAKE_CXX_FLAGS is ${CMAKE_CXX_FLAGS}")
message("CMAKE_CXX_COMPILER_ID is ${CMAKE_CXX_COMPILER_ID}")
message("CMAKE_CXX_COMPILER_VERSION is ${CMAKE_CXX_COMPILER_VERSION}")