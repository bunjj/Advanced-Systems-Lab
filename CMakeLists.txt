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

set(files 
    src/main.cpp 
    src/timing.cpp 
    src/util.cpp

    src/impl_ref/impl.cpp
    src/impl_ref/geometry.cpp
    src/impl_ref/scene.cpp

    src/impl_opt0/impl.cpp
    
    src/impl_opt1/impl.cpp
    src/impl_opt1/geometry.cpp
    src/impl_opt1/scene.cpp

    src/impl_opt3/impl.cpp
    src/impl_opt3/geometry.cpp
    src/impl_opt3/scene.cpp

    src/impl_opt4/impl.cpp
    src/impl_opt4/geometry.cpp
    src/impl_opt4/scene.cpp

    src/impl_opt5/impl.cpp
    src/impl_opt5/geometry.cpp
    src/impl_opt5/scene.cpp

    src/impl_vec1/impl.cpp
    src/impl_vec1/geometry.cpp
    src/impl_vec1/scene.cpp

    src/impl_vec2/impl.cpp
    src/impl_vec2/geometry.cpp
    src/impl_vec2/scene.cpp

    src/impl_vec3/impl.cpp
    src/impl_vec3/geometry.cpp
    src/impl_vec3/scene.cpp
    )

add_executable(${name} ${files})
target_include_directories(${name} PUBLIC include)
target_include_directories(${name} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/ext)
target_compile_options(${name} PUBLIC -Wall -Wextra -Wpedantic -mavx2 -mfma $<$<CONFIG:Debug>:-g -O0 -DDEBUG> $<$<BOOL:${SANITIZE}>:-fsanitize=undefined,address>)
target_link_options(${name} PUBLIC $<$<BOOL:${SANITIZE}>:-fsanitize=undefined,address>)
add_compile_definitions($<$<BOOL:${INSTRUMENT}>:DO_INSTRUMENT> $<$<CONFIG:Release>:DBG_MACRO_DISABLE> DBG_MACRO_NO_WARNING)
