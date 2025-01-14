if(TARGET nanothread::nanothread)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
  nanothread
  GITHUB_REPOSITORY mitsuba-renderer/nanothread
  GIT_TAG 9073b959f02da3395cdae8ed7e0e0f86b1c2ddb8
)
FetchContent_MakeAvailable(nanothread)
add_library(nanothread::nanothread ALIAS nanothread)