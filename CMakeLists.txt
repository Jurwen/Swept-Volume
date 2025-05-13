cmake_minimum_required(VERSION 3.28)

if(CMAKE_VERSION VERSION_GREATER_EQUAL "4.0.0")
    message(WARNING "CMake version is ${CMAKE_VERSION}, setting stuff for dependencies.")
    SET(CMAKE_POLICY_VERSION_MINIMUM ${REQUIRED_CMAKE_VERSION})
endif()

project(swept_volume LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)

set(FAST_ARRANGEMENT_ENABLE_PARALLEL OFF CACHE BOOL "" FORCE)
set(GEO_ENABLE_THREADS               OFF CACHE BOOL "" FORCE)

include(cli11)
include(Eigen3)
include(nlohmann-json)
Include(Marching4DLib)
include(nanothread)
include(mtet)
include(mtetcol)
include(implicit_functions)
include(convex_hull_membership)
include(unordered_dense)
include(smallvector)
include(libigl)
include(spdlog)
include(arrangement)

# project library
target_compile_options(mtetcol PRIVATE -Wno-shorten-64-to-32 -Wno-logical-op-parentheses)
set(TBB_STRICT           OFF CACHE BOOL   "" FORCE)
set(TBB_WARNING_SUPPRESS "-Wno-shorten-64-to-32" CACHE STRING "" FORCE)
file(GLOB SRC_FILES ${CMAKE_CURRENT_LIST_DIR}/src/*.cpp ${CMAKE_CURRENT_LIST_DIR}/src/*.h)
add_library(adaptive_column_grid STATIC ${SRC_FILES})
target_include_directories(adaptive_column_grid PUBLIC
     ${CMAKE_CURRENT_LIST_DIR}/src)
target_link_libraries(adaptive_column_grid
    PUBLIC 
    Eigen3::Eigen
    mtet::mtet 
    mtetcol::mtetcol
    nanothread::nanothread
    convex_hull_membership
    implicit_functions::implicit_functions
    nlohmann_json::nlohmann_json
    smallvector::smallvector
    unordered_dense::unordered_dense
    igl::core
    spdlog::spdlog
    )


add_executable(general_sweep "app/general_sweep.cpp")

target_compile_features(general_sweep PRIVATE cxx_std_17)
target_include_directories(general_sweep PUBLIC ${CMAKE_CURRENT_LIST_DIR}/src)
target_link_libraries(general_sweep PRIVATE
    adaptive_column_grid
    CLI11::CLI11
    Marching4DLib
    arrangement::arrangement
    )

option(GEN_SWEEP_TESTS "general sweep unit tests" ON)
if (GEN_SWEEP_TESTS)
 include(CTest)
 enable_testing()
 include(Catch2)

 file(GLOB TEST_FILES "${CMAKE_CURRENT_LIST_DIR}/tests/*.cpp")
 add_executable(sweep_volume_test ${TEST_FILES})
 target_link_libraries(sweep_volume_test
    PRIVATE 
    adaptive_column_grid
    Catch2::Catch2
    mtet::mtet
    mtetcol::mtetcol
    unordered_dense::unordered_dense
    smallvector::smallvector)
 target_compile_features(sweep_volume_test PRIVATE cxx_std_20)
 target_compile_definitions(sweep_volume_test PRIVATE CATCH_CONFIG_ENABLE_BENCHMARKING)
target_compile_definitions(
    sweep_volume_test
    PRIVATE                   
    -DTEST_FILE="${CMAKE_CURRENT_LIST_DIR}/data"
)
 catch_discover_tests(sweep_volume_test)
endif()
