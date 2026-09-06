# Modified in this branch, 2026.
# This build intentionally contains only the Stim-free native matching core.

cmake_minimum_required(VERSION 3.16)
project(PyMatchingCore VERSION 2.3.1 LANGUAGES CXX)

include(CTest)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(PYMATCHING_CORE_SOURCES
    src/pymatching/core/matching_region.cc
    src/pymatching/sparse_blossom/flooder/graph.cc
    src/pymatching/sparse_blossom/flooder/detector_node.cc
    src/pymatching/sparse_blossom/flooder/graph_fill_region.cc
    src/pymatching/sparse_blossom/flooder/graph_flooder.cc
    src/pymatching/sparse_blossom/flooder/match.cc
    src/pymatching/sparse_blossom/flooder_matcher_interop/compressed_edge.cc
    src/pymatching/sparse_blossom/flooder_matcher_interop/mwpm_event.cc
    src/pymatching/sparse_blossom/flooder_matcher_interop/region_edge.cc
    src/pymatching/sparse_blossom/matcher/alternating_tree.cc
    src/pymatching/sparse_blossom/matcher/mwpm.cc
    src/pymatching/sparse_blossom/search/search_detector_node.cc
    src/pymatching/sparse_blossom/search/search_flooder.cc
    src/pymatching/sparse_blossom/search/search_graph.cc
    src/pymatching/sparse_blossom/tracker/flood_check_event.cc)

add_library(libpymatching_core STATIC ${PYMATCHING_CORE_SOURCES})
set_target_properties(libpymatching_core PROPERTIES OUTPUT_NAME pymatching_core POSITION_INDEPENDENT_CODE ON)
target_include_directories(libpymatching_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
if(NOT MSVC)
    target_compile_options(libpymatching_core PRIVATE -O3 -fno-strict-aliasing)
endif()

install(TARGETS libpymatching_core ARCHIVE DESTINATION lib)
install(FILES src/pymatching/core/matching_region.h DESTINATION include/pymatching/core)

if(BUILD_TESTING)
    add_executable(pymatching_core_tests src/pymatching/core/matching_region.test.cc)
    target_link_libraries(pymatching_core_tests PRIVATE libpymatching_core)
    # Contract tests use assert, including in Release builds.
    if(MSVC)
        target_compile_options(pymatching_core_tests PRIVATE /UNDEBUG)
    else()
        target_compile_options(pymatching_core_tests PRIVATE -UNDEBUG)
    endif()
    add_test(NAME pymatching_core_tests COMMAND pymatching_core_tests)
endif()
