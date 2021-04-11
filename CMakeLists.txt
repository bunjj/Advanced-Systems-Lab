cmake_minimum_required(VERSION 3.16.0 FATAL_ERROR)
project(ASL C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
  message(STATUS "No build type specified; using ${CMAKE_BUILD_TYPE}")
endif()

set(name "main")

find_package(nlohmann_json REQUIRED)

option(INSTRUMENT "Instrument code with flop counters" OFF)

add_executable(${name} src/main.cpp src/loader.cpp src/geometry.cpp)
target_include_directories(${name} PUBLIC include)
target_compile_options(${name} PUBLIC -Wall -Wextra -Wpedantic $<$<CONFIG:Debug>:-fsanitize=undefined,address -g -O0 -DDEBUG>)
target_link_options(${name} PUBLIC $<$<CONFIG:Debug>:-fsanitize=undefined,address>)
target_link_libraries(${name} PRIVATE nlohmann_json::nlohmann_json)

if (INSTRUMENT)
    add_compile_definitions(DO_INSTRUMENT) 
endif()
