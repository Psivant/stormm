#--- Define the test application targets ---

set(STORMM_TEST_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/UnitTesting/test_unit_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/FileManagement/test_file_management.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Accelerator/test_hybrid.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Parsing/test_parse.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Parsing/test_input.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Parsing/test_output.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Namelists/test_namelists.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Math/test_clustering.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Math/test_math.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Math/test_numerics.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Topology/test_amber_prmtop.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Topology/test_topology_mods.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Trajectory/test_amber_coordinates.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Trajectory/test_coordinate_swap.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Trajectory/test_coordinate_formats.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Chemistry/test_chemical_features.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Chemistry/test_atommask.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Structure/test_local_arrangement.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Structure/test_isomerization.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Structure/test_mesh_construction.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Sampling/test_enhanced_sampling.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Synthesis/test_atomgraph_synthesis.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Synthesis/test_condensate.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Synthesis/test_small_molecules.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Synthesis/test_synthesis.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Potential/test_valence_evaluation.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Potential/test_nonbonded_evaluation.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Potential/test_generalized_born.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Potential/test_neighbor_list.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Potential/test_pme.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Reporting/test_display.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Reporting/test_progressbar.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Restraints/test_restraints.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/MolecularMechanics/test_dynamics.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/MolecularMechanics/test_minimization.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/MoleculeFormat/test_molecule_formats.cpp)

if (STORMM_ENABLE_RDKIT)

    list(APPEND STORMM_TEST_SOURCES
            ${CMAKE_CURRENT_SOURCE_DIR}/RDKit/test_rdkit.cpp)

endif()

if(STORMM_ENABLE_CUDA)

    list(APPEND STORMM_TEST_SOURCES
            ${CMAKE_CURRENT_SOURCE_DIR}/Accelerator/test_hpc_hybrid.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/Accelerator/test_hpc_status.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/Math/test_hpc_math.cu
            ${CMAKE_CURRENT_SOURCE_DIR}/MolecularMechanics/test_hpc_minimization.cu
            ${CMAKE_CURRENT_SOURCE_DIR}/MolecularMechanics/test_hpc_dynamics.cu
            ${CMAKE_CURRENT_SOURCE_DIR}/Synthesis/test_hpc_synthesis.cu)

endif()

foreach (STORMM_TEST_SOURCE ${STORMM_TEST_SOURCES})

    get_filename_component(STORMM_TEST_NAME ${STORMM_TEST_SOURCE} NAME_WE)

    add_executable(${STORMM_TEST_NAME} ${STORMM_TEST_SOURCE})
    target_link_libraries(${STORMM_TEST_NAME} ${PROJECT_NAME})

    if(STORMM_TEST_NAME STREQUAL "test_atomgraph_synthesis" OR
       STORMM_TEST_NAME STREQUAL "test_hpc_synthesis")
        # Add extra dependency for the `test_atomgraph_synthesis` test.
        target_sources(${STORMM_TEST_NAME} PRIVATE
                ${CMAKE_CURRENT_SOURCE_DIR}/Synthesis/assemble_restraints.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/Synthesis/assemble_restraints.h)
    endif()

    if (STORMM_ENABLE_RDKIT)
        target_link_libraries(${STORMM_TEST_NAME} RDKit::SmilesParse)
    endif()
    if(STORMM_ENABLE_CUDA)
        target_include_directories(${STORMM_TEST_NAME} PRIVATE ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
        target_link_libraries(${STORMM_TEST_NAME} ${STORMM_CUDA_LIBS})
    endif()

    # Register the test with CTest (and hence also `make test`)
    add_test(NAME ${STORMM_TEST_NAME} COMMAND ${STORMM_TEST_NAME} WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})

endforeach()
