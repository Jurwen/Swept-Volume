if (TARGET mtetcol::mtetcol)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
    mtetcol
    GIT_REPOSITORY git@github.com:qnzhou/mtetcol.git
    GIT_TAG v0.0.6
    )

FetchContent_MakeAvailable(mtetcol)
