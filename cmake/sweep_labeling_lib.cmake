if(TARGET sweep_labeling_lib)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
    sweep_labeling_lib
    GIT_REPOSITORY https://github.com/duxingyi-charles/SweepLabeling.git
    GIT_TAG f975afbf3b1261ca346c5fb5fb479d1a6e1505b9
)
FetchContent_MakeAvailable(sweep_labeling_lib)
