cmake_minimum_required(VERSION 3.10.1)

#=== Set name of our project
project(project3)

set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_EXTENSIONS OFF)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)

#set(CMAKE_BUILD_TYPE Debug)

#set(CMAKE_CXX_FLAGS "-Wall -Werror -pedantic -Wextra -g -Og -Oz -fsanitize=address")
#set(CMAKE_CXX_FLAGS "-Wall -Werror -Wextra -Wpedantic -O0 -g3 -Og -glldb -Oz -fsanitize=address -fverbose-asm")
set(CMAKE_CXX_FLAGS "-Wall -Werror -Wextra -Wpedantic -g3")
#set(CMAKE_CXX_FLAGS "-Wall -Werror -pedantic -Wextra -Ofast -O3")
#set(CMAKE_CXX_FLAGS "-Wall -Werror -pedantic -Wextra -g -fverbose-asm")
#set(CMAKE_CXX_FLAGS "-Wall -Wextra -g -fverbose-asm -masm=intel")
#set(CMAKE_CXX_FLAGS "-Wall -Werror -pedantic -Wextra -O3")
#set(CMAKE_CXX_FLAGS "-Wall -Wextra -g -c -Wa,-alh -masm=intel -fverbose-asm")

# Always produce a shared library if we do not pass a library type to the add_library command
#set(BUILD_SHARED_LIBS ON)

#find_package(OpenMP REQUIRED)
# Bring headers into the project
include_directories(include)
include_directories(include/Ar)

#=== Find all source files in a directory
aux_source_directory(./src SOURCE_FILES)

# Add Sources
file( GLOB SOURCES "src/*")

set( LIBRARY_OUTPUT_PATH "${CMAKE_CURRENT_SOURCE_DIR}/build/lib" )
add_library( project3 SHARED ${SOURCES})
#add_library(project3 ${SOURCES})

set( EXECUTABLE_OUTPUT_PATH "${CMAKE_CURRENT_SOURCE_DIR}/build/scripts")
add_executable( main scripts/main.cxx )
target_link_libraries(main project3)

#message("CMAKE_BUILD_TYPE is ${CMAKE_BUILD_TYPE}")
message("CMAKE_CXX_FLAGS is ${CMAKE_CXX_FLAGS}")

#message("CMAKE_C_FLAGS_DEBUG is ${CMAKE_C_FLAGS_DEBUG}")
#message("CMAKE_C_FLAGS_RELEASE is ${CMAKE_C_FLAGS_RELEASE}")
#message("CMAKE_C_FLAGS_RELWITHDEBINFO is ${CMAKE_C_FLAGS_RELWITHDEBINFO}")
#message("CMAKE_C_FLAGS_MINSIZEREL is ${CMAKE_C_FLAGS_MINSIZEREL}")

#message("CMAKE_CXX_FLAGS_DEBUG is ${CMAKE_CXX_FLAGS_DEBUG}")
#message("CMAKE_CXX_FLAGS_RELEASE is ${CMAKE_CXX_FLAGS_RELEASE}")
#message("CMAKE_CXX_FLAGS_RELWITHDEBINFO is ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
#message("CMAKE_CXX_FLAGS_MINSIZEREL is ${CMAKE_CXX_FLAGS_MINSIZEREL}")

message("CMAKE_CXX_COMPILER_ID is ${CMAKE_CXX_COMPILER_ID}")
message("CMAKE_CXX_COMPILER_VERSION is ${CMAKE_CXX_COMPILER_VERSION}")