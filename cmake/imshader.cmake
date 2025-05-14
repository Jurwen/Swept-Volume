if (TARGET implicit_shader::imshader)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
    imshader
    GIT_REPOSITORY https://github.com/qnzhou/implicit_shader.git
    GIT_TAG main
    )

FetchContent_MakeAvailable(imshader)
