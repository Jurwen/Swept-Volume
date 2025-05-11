if(TARGET igl::core)
    return()
endif()

message(STATUS "Third-party: creating target 'igl::core'")
set(LIBIGL_COPYLEFT_CGAL ON  CACHE BOOL "" FORCE)
set(LIBIGL_INSTALL       OFF CACHE BOOL "" FORCE)
set(Boost_NO_BOOST_CMAKE  ON   CACHE BOOL "" FORCE)
set(BOOST_ROOT            "/opt/homebrew/opt/boost" CACHE PATH "" FORCE)
include(FetchContent)
CPMAddPackage(
    libigl
    GIT_REPOSITORY https://github.com/libigl/libigl.git
    GIT_TAG 3ea7f9480967fcf6bf02ce9b993c0ea6d2fc45f6
    OPTIONS 
        LIBIGL_COPYLEFT_CGAL:BOOL=ON
        LIBIGL_INSTALL OFF
)
# include(eigen)
FetchContent_MakeAvailable(libigl)
