cmake_minimum_required(VERSION 2.6)
project(splitter)

link_directories(~/cpp/htslib/)
include_directories(~/cpp/htslib/)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

file(GLOB splitter_SRC
    "GffLoader/*.h"
    "GffLoader/*.cpp"
    "SplitFinder/*.h"
    "SplitFinder/*.cpp"
    "Utils/*.h"
    "Utils/*.cpp"
    "src/*.h"
    "src/*.cpp"
)

add_executable(splitter ${splitter_SRC})

target_link_libraries (splitter hts)

install(TARGETS splitter RUNTIME DESTINATION bin)
