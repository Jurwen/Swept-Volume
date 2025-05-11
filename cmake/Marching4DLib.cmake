if (TARGET Marching4DLib)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
    Marching4DLib
    GIT_REPOSITORY git@github.com:duxingyi-charles/Marching4D.git
    GIT_TAG 11a48c25a1041649bc152ea171ecc27b4937abb3
    )

FetchContent_MakeAvailable(Marching4DLib)
