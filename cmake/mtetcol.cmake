if (TARGET mtetcol::mtetcol)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
    mtetcol
    GIT_REPOSITORY git@github.com:qnzhou/mtetcol.git
    GIT_TAG main
    )

FetchContent_MakeAvailable(mtetcol)
