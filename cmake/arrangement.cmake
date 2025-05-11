if(TARGET arrangement::arrangement)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
    arrangement
    GIT_REPOSITORY https://github.com/qnzhou/arrangement-benchmark.git
    GIT_TAG 28fbb6abcd24401861376020c53973998de30089
)
FetchContent_MakeAvailable(arrangement)
