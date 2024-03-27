cmake_minimum_required(VERSION 3.18)

if (CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    message(FATAL_ERROR "CMake should not be called from the root directory.")
endif()

project(stormm LANGUAGES CXX)
project(stormm VERSION 0.1)
string(TOUPPER ${PROJECT_NAME} PROJECT_NAME_UPPER)

#--- Define any options the user can set ---

option(STORMM_BUILD_APPS "Build the apps" ON)
option(STORMM_BUILD_TESTS "Build the tests" ON)
option(STORMM_BUILD_BENCHMARKS "Build the benchmarks" ON)

option(STORMM_ENABLE_CUDA "Enable CUDA support" ON)
option(STORMM_ENABLE_RDKIT "Enable RDKIT support" ON)

option(STORMM_ENABLE_TEST_COVERAGE "Enable test coverage" OFF)

#--- Define any global settings / flags ---

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)  # set the -fPIC flag

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -O0 -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3")
set(STORMM_VERSION "0.3")

if (STORMM_ENABLE_RDKIT) 
    message("-- RDKit support will be enabled")

    find_package(RDKit REQUIRED)

    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -DSTORMM_USE_RDKIT -DSTORMM_USE_RDKIT")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSTORMM_USE_RDKIT -DSTORMM_USE_RDKIT")
endif()

if (STORMM_ENABLE_CUDA)
    message("-- CUDA support will be enabled")

    enable_language(CUDA)

    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)

    set(CMAKE_CUDA_FLAGS_DEBUG "${CMAKE_CUDA_FLAGS_DEBUG} --compiler-options=-O0 -g")

    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -DSTORMM_USE_HPC -DSTORMM_USE_CUDA --ptxas-options=\"-v\"")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSTORMM_USE_HPC -DSTORMM_USE_CUDA")

    if (NOT DEFINED CUSTOM_GPU_ARCH)
        set(CMAKE_CUDA_ARCHITECTURES "52;60;61;70;75;80;86;89")
    else()
        set(CMAKE_CUDA_ARCHITECTURES "${CUSTOM_GPU_ARCH}")
    endif()

    find_package(CUDAToolkit REQUIRED)
    set(STORMM_CUDA_LIBS CUDA::curand CUDA::cublas CUDA::cusolver CUDA::cudart CUDA::nvml CUDA::cufft)

endif()

if (STORMM_ENABLE_TEST_COVERAGE)
    if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "(Apple)?[Cc]lang")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-instr-generate -fcoverage-mapping")
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage -fprofile-arcs -ftest-coverage")
    endif()
endif()

# --- Create a version file ---

## We pull the version information directly from Git tags, so that automating releases is much
## easier.

find_package(Git)

if(GIT_EXECUTABLE)
    execute_process(
        COMMAND ${GIT_EXECUTABLE} describe --tags --dirty
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_DESCRIBE_VERSION
        RESULT_VARIABLE GIT_DESCRIBE_ERROR_CODE
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(NOT GIT_DESCRIBE_ERROR_CODE)
        set(STORMM_VERSION ${GIT_DESCRIBE_VERSION})
    endif()
endif()

if(NOT DEFINED STORMM_VERSION)
    set(STORMM_VERSION 0+unknown)
    message(WARNING "Could not find STORMM_VERSION using git... using default version \"${STORMM_VERSION}\".")
endif()

configure_file(
        "${CMAKE_SOURCE_DIR}/src/version.h.in" 
        "${CMAKE_BINARY_DIR}/version.h" 
        @ONLY)

# --- Define any source files ---

## Using GLOB is not considered a 'best practice' as new files aren't always detected by CMake's
## incremental build system, and it may also include any extra source files in the tree that you
## may not want to have included... but is still a very widely used pattern to avoid the need to
## manually list all files associated with the library making it much easier to use and maintain.
file(GLOB_RECURSE STORMM_CPP_FILES CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/src/*.cpp")
file(GLOB_RECURSE STORMM_TPP_FILES CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/src/*.tpp")
file(GLOB_RECURSE STORMM_HPP_FILES CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/src/*.h")

file(GLOB_RECURSE STORMM_CUDA_FILES          CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/src/*.cu")
file(GLOB_RECURSE STORMM_CUDA_HEADERS        CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/src/*.cuh")
file(GLOB_RECURSE STORMM_CUDA_INCLUDED_FILES CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/src/*.cui")

# --- Define the library ---

add_library(${PROJECT_NAME} SHARED ${STORMM_CPP_FILES} ${STORMM_TPP_FILES})

target_include_directories(${PROJECT_NAME}
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
        $<INSTALL_INTERFACE:include/${PROJECT_NAME}>)

if (STORMM_ENABLE_CUDA)
    target_sources(${PROJECT_NAME} PRIVATE ${STORMM_CUDA_FILES} ${STORMM_CUDA_HEADERS} ${STORMM_CUDA_INCLUDED_FILES})
    target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
    target_link_libraries(${PROJECT_NAME} ${STORMM_CUDA_LIBS})
endif()

if (STORMM_ENABLE_RDKIT)
    target_link_libraries(${PROJECT_NAME} RDKit::SmilesParse)
endif()

# --- Create an installable target ---

include(GNUInstallDirs)   # Define 'sensible' defaults for installation directories on UNIX

install(DIRECTORY "${CMAKE_SOURCE_DIR}/src/"
        DESTINATION "include/${PROJECT_NAME}"
        FILES_MATCHING PATTERN "*.h")

install(TARGETS ${PROJECT_NAME}
        EXPORT  ${PROJECT_NAME_UPPER}Targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        BUNDLE  DESTINATION ${CMAKE_INSTALL_BINDIR})

install(EXPORT ${PROJECT_NAME_UPPER}Targets
        FILE   ${PROJECT_NAME_UPPER}Targets.cmake
        NAMESPACE STORMM::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})        

#--- Build the tests ---
if (STORMM_BUILD_TESTS)
    include(CTest)
    add_subdirectory(test)
endif()

#--- Build the apps ---
if (STORMM_BUILD_APPS)
    add_subdirectory(apps)
endif()

#--- Build the benchmarks ---
if (STORMM_BUILD_BENCHMARKS)
    add_subdirectory(benchmark)
endif()  
