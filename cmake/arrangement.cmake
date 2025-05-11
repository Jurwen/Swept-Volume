if(TARGET arrangement::arrangement)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
    arrangement
    GIT_REPOSITORY https://github.com/qnzhou/arrangement-benchmark.git
    GIT_TAG main
)
FetchContent_MakeAvailable(arrangement)
