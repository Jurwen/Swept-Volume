if (TARGET Marching4DLib)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
    Marching4DLib
    GIT_REPOSITORY https://github.com/duxingyi-charles/Marching4D.git
    GIT_TAG main
    )

FetchContent_MakeAvailable(Marching4DLib)
