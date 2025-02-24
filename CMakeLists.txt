cmake_minimum_required(VERSION 3.15)

project(swept_volume LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)

include(cli11)
include(Eigen3)
include(nlohmann-json)
Include(Marching4DLib)
include(nanothread)
include(mtet)
include(implicit_functions)
include(convex_hull_membership)
include(unordered_dense)
include(smallvector)

# libigl
option(LIBIGL_WITH_CGAL "Use CGAL"          ON)
find_package(LIBIGL REQUIRED QUIET)


# project library
file(GLOB SRC_FILES ${CMAKE_CURRENT_LIST_DIR}/src/*.cpp ${CMAKE_CURRENT_LIST_DIR}/src/*.h)
add_library(adaptive_column_grid STATIC ${SRC_FILES})
target_include_directories(adaptive_column_grid PUBLIC
     ${CMAKE_CURRENT_LIST_DIR}/src)
target_link_libraries(adaptive_column_grid
    PUBLIC 
    Eigen3::Eigen
    mtet::mtet 
    nanothread::nanothread
    convex_hull_membership
    implicit_functions::implicit_functions
    nlohmann_json::nlohmann_json
    smallvector::smallvector
    unordered_dense::unordered_dense
    igl::core
    igl::cgal
    )

add_executable(general_sweep "app/general_sweep")

target_compile_features(general_sweep PRIVATE cxx_std_17)
target_include_directories(general_sweep PUBLIC ${CMAKE_CURRENT_LIST_DIR}/src)
target_link_libraries(general_sweep PRIVATE
    adaptive_column_grid
    CLI11::CLI11
    Marching4DLib
    )
