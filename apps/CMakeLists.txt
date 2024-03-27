#--- Build the apps ---

set(STORMM_APP_SUFFIX "stormm")

if (STORMM_ENABLE_CUDA)
    set(STORMM_APP_SUFFIX "${STORMM_APP_SUFFIX}.cuda")
endif()

add_subdirectory(Dyna)
add_subdirectory(Ffrn)
add_subdirectory(Conf)
