cmake_minimum_required(VERSION 3.16.0 FATAL_ERROR)
project(ASL C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
  message_colored(STATUS "No build type specified; using ${CMAKE_BUILD_TYPE}" 33)
endif()

set(name "main")

find_package(nlohmann_json REQUIRED)

add_executable(${name} src/main.cpp src/loader.cpp)
target_include_directories(${name} PUBLIC include)
target_compile_options(${name} PUBLIC -Wall -Wextra -Wpedantic $<$<CONFIG:Instrument>:-fsanitize=undefined,address -g -O0 -DDEBUG>)
target_link_options(${name} PUBLIC $<$<CONFIG:Instrument>:-fsanitize=undefined,address>)
target_link_libraries(${name} PRIVATE nlohmann_json::nlohmann_json)
add_compile_definitions($<$<CONFIG:Instrument>:DO_INSTRUMENT>)
