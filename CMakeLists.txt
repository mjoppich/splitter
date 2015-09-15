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
    "Utils/*.h"
    "Utils/*.cpp"
    "GffLoader/*.h"
    "GffLoader/*.cpp"
    "SplitFinder/*.h"
    "SplitFinder/*.cpp"
)

file(GLOB gtftools_SRC
    "Utils/*.h"
    "Utils/*.cpp"
    "GffLoader/*.h"
    "GffLoader/*.cpp"
)

set(CMAKE_CXX_STANDARD 11)

add_executable(splitter ${splitter_SRC} "src/main.cpp")
add_executable(gtftool ${gtftools_SRC} "src/mainGTFTool.cpp")

target_link_libraries (splitter hts)

install(TARGETS splitter gtftool RUNTIME DESTINATION bin)
