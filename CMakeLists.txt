cmake_minimum_required(VERSION 3.28)

option(SWEEP_PYTHON_BINDING "Build python binding" OFF)

if(CMAKE_VERSION VERSION_GREATER_EQUAL "4.0.0")
    message(WARNING "CMake version is ${CMAKE_VERSION}, setting stuff for dependencies.")
    SET(CMAKE_POLICY_VERSION_MINIMUM ${REQUIRED_CMAKE_VERSION})
endif()

project(swept_volume LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POLICY_VERSION_MINIMUM 3.28)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)

set(FAST_ARRANGEMENT_ENABLE_PARALLEL OFF CACHE BOOL "" FORCE)
set(GEO_ENABLE_THREADS               OFF CACHE BOOL "" FORCE)

include(cli11)
include(Eigen3)
include(nlohmann-json)
include(mtet)
include(mtetcol)
include(mshio)
include(unordered_dense)
include(smallvector)
include(libigl)
include(spdlog)
include(arrangement)
include(stf)
include(lagrange)
include(yaml-cpp)
include(nanothread)

# project library
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|AppleClang|GNU")
    target_compile_options(mtetcol PRIVATE -Wno-shorten-64-to-32 -Wno-logical-op-parentheses)
endif()
set(TBB_STRICT           OFF CACHE BOOL   "" FORCE)
set(TBB_WARNING_SUPPRESS "-Wno-shorten-64-to-32" CACHE STRING "" FORCE)
file(GLOB SRC_FILES ${CMAKE_CURRENT_LIST_DIR}/src/*.cpp ${CMAKE_CURRENT_LIST_DIR}/src/*.h)
add_library(sweep STATIC ${SRC_FILES})
target_include_directories(sweep
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/include
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/src
)
target_link_libraries(sweep
    PUBLIC 
        Eigen3::Eigen
        lagrange::core
    PRIVATE
        arrangement::arrangement
        mshio::mshio
        mtet::mtet 
        mtetcol::mtetcol
        nlohmann_json::nlohmann_json
        smallvector::smallvector
        spdlog::spdlog
        unordered_dense::unordered_dense
        nanothread::nanothread
        yaml-cpp::yaml-cpp
        stf::stf
)
add_library(sweep::sweep ALIAS sweep)
set_target_properties(sweep PROPERTIES POSITION_INDEPENDENT_CODE ON)

# Optional: enable geometric verification of envelope_edge_id inheritance.
# Adds a second pass in compute_envelope_arrangement that checks each valid
# label geometrically and probes missed inheritances. Intended for testing
# during development of the edge-label feature; leave OFF for normal builds.
option(SWEEP_VERIFY_EDGE_LABELS "Enable geometric verification of edge-label inheritance" OFF)
if(SWEEP_VERIFY_EDGE_LABELS)
    target_compile_definitions(sweep PRIVATE SWEEP_VERIFY_EDGE_LABELS)
endif()

file(GLOB SWEEP_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/app/include/*.h)
add_executable(generalized_sweep "app/generalized_sweep.cpp")

target_compile_features(generalized_sweep PRIVATE cxx_std_17)
target_compile_definitions(generalized_sweep PRIVATE -DDATA_DIR="${CMAKE_CURRENT_LIST_DIR}/data")
target_include_directories(generalized_sweep PRIVATE ${CMAKE_CURRENT_LIST_DIR}/app/include)
target_sources(generalized_sweep PRIVATE ${SWEEP_HEADERS})
set_target_properties(generalized_sweep PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(generalized_sweep
    PRIVATE
        sweep::sweep
        CLI11::CLI11
        stf::stf
        lagrange::core
        lagrange::io
        nlohmann_json::nlohmann_json
        yaml-cpp::yaml-cpp
)

# On Windows there is no RPATH, so shared-library dependencies (e.g. nanothread.dll)
# must sit next to the executable at runtime. Copy all runtime DLLs alongside the
# exe after every build. No-op on Linux/macOS, where it evaluates to an empty list.
if(WIN32)
    add_custom_command(TARGET generalized_sweep POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_RUNTIME_DLLS:generalized_sweep>
            $<TARGET_FILE_DIR:generalized_sweep>
        COMMAND_EXPAND_LISTS
        VERBATIM
    )
endif()

option(GEN_SWEEP_TESTS "general sweep unit tests" OFF)
if (GEN_SWEEP_TESTS)
 include(CTest)
 enable_testing()
 include(Catch2)

 file(GLOB TEST_FILES "${CMAKE_CURRENT_LIST_DIR}/tests/*.cpp")
 add_executable(sweep_volume_test ${TEST_FILES})
 target_link_libraries(sweep_volume_test
    PRIVATE 
    sweep
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

 if(WIN32)
     add_custom_command(TARGET sweep_volume_test POST_BUILD
         COMMAND ${CMAKE_COMMAND} -E copy_if_different
             $<TARGET_RUNTIME_DLLS:sweep_volume_test>
             $<TARGET_FILE_DIR:sweep_volume_test>
         COMMAND_EXPAND_LISTS
         VERBATIM
     )
 endif()
endif()

if (SWEEP_PYTHON_BINDING)
    if(SKBUILD)
        if(APPLE)
            set(CMAKE_INSTALL_RPATH @loader_path)
        elseif(UNIX)
            set(CMAKE_INSTALL_RPATH $ORIGIN)
        endif()
    endif()

    include(nanobind)
    set(PY_SRC_FILE ${CMAKE_CURRENT_LIST_DIR}/python/pysweep3d.cpp)
    nanobind_add_module(pysweep3d NB_STATIC ${PY_SRC_FILE})
    target_link_libraries(pysweep3d PRIVATE sweep::sweep nanothread::nanothread)
    target_compile_features(pysweep3d PRIVATE cxx_std_20)
    install(TARGETS pysweep3d nanothread LIBRARY
        COMPONENT Sweep_Python_Runtime
        DESTINATION ${SKBUILD_PLATLIB_DIR}/sweep3d
    )
    add_library(sweep::pysweep3d ALIAS pysweep3d)
endif()
