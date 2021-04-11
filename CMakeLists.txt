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

option(INSTRUMENT "Instrument code with flop counters" OFF)
option(SANITIZE "Adds UB and Address Sanitizers" OFF)

add_executable(${name} src/main.cpp src/loader.cpp src/geometry.cpp)
target_include_directories(${name} PUBLIC include)
target_include_directories(${name} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/ext)
target_compile_options(${name} PUBLIC -Wall -Wextra -Wpedantic $<$<CONFIG:Debug>:-g -O0 -DDEBUG> $<$<BOOL:${SANITIZE}>:-fsanitize=undefined,address>)
target_link_options(${name} PUBLIC $<$<BOOL:${SANITIZE}>:-fsanitize=undefined,address>)
add_compile_definitions($<$<BOOL:${INSTRUMENT}>:DO_INSTRUMENT> $<$<CONFIG:Release>:DBG_MACRO_DISABLE> DBG_MACRO_NO_WARNING)
