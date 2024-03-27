#--- Define the test application targets ---

set(STORMM_BENCHMARK_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/DirectSpace/hail_examination.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/DirectSpace/erfc_tables.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/ReciprocalSpace/charge_mapping.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/ForceAccumulation/valence.cpp)

if(STORMM_ENABLE_CUDA)

    list(APPEND STORMM_BENCHMARK_SOURCES
            ${CMAKE_CURRENT_SOURCE_DIR}/ForceAccumulation/accumulate.cu
            ${CMAKE_CURRENT_SOURCE_DIR}/KernelTesting/test_nonperiodic_kernels.cu
            ${CMAKE_CURRENT_SOURCE_DIR}/ReciprocalSpace/fft_throughput.cu)

endif()

foreach (STORMM_BENCHMARK_SOURCE ${STORMM_BENCHMARK_SOURCES})

    get_filename_component(STORMM_BENCHMARK_NAME ${STORMM_BENCHMARK_SOURCE} NAME_WE)

    add_executable(${STORMM_BENCHMARK_NAME} ${STORMM_BENCHMARK_SOURCE})
    target_link_libraries(${STORMM_BENCHMARK_NAME} ${PROJECT_NAME})

    if (STORMM_ENABLE_RDKIT)
        target_link_libraries(${STORMM_BENCHMARK_NAME} RDKit::SmilesParse)
    endif()
    if(STORMM_ENABLE_CUDA)
        target_include_directories(${STORMM_BENCHMARK_NAME} PRIVATE ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
        target_link_libraries(${STORMM_BENCHMARK_NAME} ${STORMM_CUDA_LIBS})
    endif()

    # Register the test with CTest (and hence also `make benchmark`)
    #add_test(NAME ${STORMM_BENCHMARK_NAME} COMMAND ${STORMM_BENCHMARK_NAME} WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})

endforeach()
