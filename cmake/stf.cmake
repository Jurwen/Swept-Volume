if (TARGET stf::stf)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
    stf
    GIT_REPOSITORY git@github.com:qnzhou/space-time-functions.git
    GIT_TAG main
    )

FetchContent_MakeAvailable(stf)
