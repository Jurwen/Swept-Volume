if (TARGET stf::stf)
    return()
endif()

include(CPM)
CPMAddPackage(
    NAME stf
    GITHUB_REPOSITORY duxingyi-charles/space-time-functions
    GIT_TAG 46580b8
)
