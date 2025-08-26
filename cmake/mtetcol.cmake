if (TARGET mtetcol::mtetcol)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
    mtetcol
    GIT_REPOSITORY git@github.com:qnzhou/mtetcol.git
    GIT_TAG ba7951bf217f535ee1fccb93c348194e04155b64
    )

FetchContent_MakeAvailable(mtetcol)
