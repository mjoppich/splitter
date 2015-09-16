#MakeLists files in this project can
# refer to the root source directory of the project as ${HELLO_SOURCE_DIR} and
# to the root binary directory of the project as ${HELLO_BINARY_DIR}.
cmake_minimum_required (VERSION 2.8.11)
project (splittools)
link_directories(~/cpp/htslib/)
include_directories(~/cpp/htslib/)
find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
set(CMAKE_CXX_STANDARD 11)
file(GLOB splitutilslib_SRC
        "utils/*.h"
        "utils/*.cpp"
)
add_library (splitutils ${splitutilslib_SRC})
file(GLOB gtxlib_SRC
        "gtxloader/*.h"
        "gtxloader/*.cpp"
)
add_library (gtx ${gtxlib_SRC})
file(GLOB splitslib_SRC
        "splits/*.h"
        "splits/*.cpp"
)
add_library (splits ${splitslib_SRC})
TARGET_LINK_LIBRARIES(gtx splitutils)
TARGET_LINK_LIBRARIES(splits gtx hts)
ADD_EXECUTABLE(gtxtools "tools/gtxtools/main.cpp")
ADD_EXECUTABLE(splitter "tools/splitter/main.cpp")
TARGET_LINK_LIBRARIES(gtxtools gtx)
TARGET_LINK_LIBRARIES(splitter splits)
target_include_directories (splitutils PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories (gtx PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories (splits PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories (gtxtools PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories (splitter PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
