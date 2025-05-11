if(TARGET igl::core)
    return()
endif()

message(STATUS "Third-party: creating target 'igl::core'")

set(LIBIGL_COPYLEFT_CGAL On CACHE BOOL "Build CGAL support in libigl")

include(FetchContent)
CPMAddPackage(
    libigl
    GIT_REPOSITORY https://github.com/libigl/libigl.git
    GIT_TAG ba69acc509632c56c738f03440841cc37da50497
    OPTIONS LIBIGL_INSTALL OFF
)
# include(eigen)
FetchContent_MakeAvailable(libigl)
