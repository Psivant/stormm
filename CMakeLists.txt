# STORMM Build System Configuration
#
# This CMakeLists.txt orchestrates the build process for STORMM, handling both CPU and GPU components.
# The build is structured in phases to optimize memory usage:
# 1. CPU components are compiled first as object files
# 2. After CPU compilation completes, memory is freed before GPU compilation begins
# 3. Finally, all components are linked into the final library
#
# A simple 'make -j' will automatically handle these phases in the correct order.
# The build system will display messages when each phase completes, indicating when
# memory can be freed.

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

# Preserve compatibility with older build scripts while migrating to STORMM_ENABLE_PYTHON.
if (DEFINED STORMM_BUILD_PYTHON AND NOT DEFINED STORMM_ENABLE_PYTHON)
    set(STORMM_ENABLE_PYTHON ${STORMM_BUILD_PYTHON})
    message(WARNING "STORMM_BUILD_PYTHON is deprecated; use STORMM_ENABLE_PYTHON instead.")
endif()
option(STORMM_ENABLE_PYTHON "Enable the Python bridge library" ON)
option(STORMM_ENABLE_CUDA "Enable CUDA support" ON)
option(STORMM_ENABLE_RDKIT "Enable RDKIT support" ON)
option(STORMM_ENABLE_TEST_COVERAGE "Enable test coverage" OFF)
option(STORMM_INCLUDE_POCKETFFT "Enable PocketFFT support" ON)
option(STORMM_INCLUDE_NETCDF "Enable NetCDF Support" OFF)
option(STORMM_IGNORE_RAM_LIMITS "Ignore RAM limits" ON)

# --- PocketFFT Integration ---
include(ExternalProject)
if (STORMM_INCLUDE_POCKETFFT)
    message("-- Including PocketFFT as part of the build process")

    # Define the PocketFFT INTERFACE library target
    add_library(PocketFFT INTERFACE)

    # Set include directory for PocketFFT
    target_include_directories(PocketFFT INTERFACE
        ${CMAKE_BINARY_DIR}/third_party/pocketfft
    )

    # Fetch the source if not already downloaded
    set(PocketFFT_SOURCE_DIR ${CMAKE_BINARY_DIR}/third_party/pocketfft)
    if (NOT EXISTS ${PocketFFT_SOURCE_DIR})
        message("-- Cloning PocketFFT repository")
        execute_process(
            COMMAND git clone --branch cpp --depth 1 https://gitlab.mpcdf.mpg.de/mtr/pocketfft.git ${PocketFFT_SOURCE_DIR}
            RESULT_VARIABLE GIT_CLONE_RESULT
        )
        if (NOT GIT_CLONE_RESULT EQUAL 0)
            message(FATAL_ERROR "Failed to clone PocketFFT repository.")
        endif()
    endif()

    # Create a dummy build target for PocketFFT
    add_custom_target(PocketFFTBuild ALL
        COMMENT "PocketFFT is a header-only library; no build steps required."
    )

    # Add dependency between PocketFFTBuild and the PocketFFT library
    add_dependencies(PocketFFTBuild PocketFFT)

    # Link math library for users of PocketFFT
    target_link_libraries(PocketFFT INTERFACE m)

    # Install rules
    install(DIRECTORY ${PocketFFT_SOURCE_DIR}/
        DESTINATION include/pocketfft
        FILES_MATCHING PATTERN "*.h"
    )
endif()

# --- NetCDF + HDF5 Integration ---
if (STORMM_INCLUDE_NETCDF)
    message("-- Including HDF5 and NetCDF as part of the build process")

    # Set installation directories
    set(THIRD_PARTY_INSTALL_DIR ${CMAKE_BINARY_DIR}/third_party)
    set(HDF5_INSTALL_DIR ${THIRD_PARTY_INSTALL_DIR}/hdf5/install)
    set(NETCDF_INSTALL_DIR ${THIRD_PARTY_INSTALL_DIR}/netcdf)
    set(NETCDF_C_INSTALL_DIR ${NETCDF_INSTALL_DIR}/install)
    set(NETCDF_CXX_INSTALL_DIR ${NETCDF_INSTALL_DIR}/install-cxx)

    # --- Install HDF5 ---
    # --- Install HDF5 ---
    ExternalProject_Add(
      HDF5
      GIT_REPOSITORY https://github.com/HDFGroup/hdf5.git
      GIT_TAG hdf5-1_14_3
      GIT_PROGRESS 1
      GIT_SHALLOW 1
      CMAKE_ARGS
          -DCMAKE_INSTALL_PREFIX=${HDF5_INSTALL_DIR}
          -DBUILD_SHARED_LIBS=ON
          -DHDF5_BUILD_HL_LIB=ON
          -DHDF5_ENABLE_Z_LIB_SUPPORT=ON
          -DHDF5_BUILD_CPP_LIB=ON
          -DHDF5_BUILD_TOOLS=OFF
          -DHDF5_BUILD_EXAMPLES=OFF
          -DHDF5_BUILD_JAVA=OFF
          -DHDF5_BUILD_FORTRAN=OFF
          -DHDF5_BUILD_UTILS=OFF
          -DHDF5_ENABLE_PARALLEL=OFF
          -DHDF5_ENABLE_THREADSAFE=OFF
          -DHDF5_ENABLE_DIRECT_VFD=ON
          -DHDF5_ENABLE_EMBEDDED_LIBINFO=ON
          # Add these lines to disable network features that require CURL
          -DHDF5_ENABLE_ROS3_VFD=OFF
          -DHDF5_ENABLE_DAP=OFF
      BUILD_COMMAND ${CMAKE_COMMAND} --build . --parallel
      INSTALL_COMMAND ${CMAKE_COMMAND} --install .
    )

    # --- Install NetCDF-C (Requires HDF5) ---
    ExternalProject_Add(
        NetCDFC
        GIT_REPOSITORY https://github.com/Unidata/netcdf-c.git
        GIT_TAG v4.9.2
        GIT_PROGRESS 1
        DEPENDS HDF5
        CMAKE_ARGS
            -DCMAKE_INSTALL_PREFIX=${NETCDF_C_INSTALL_DIR}
            -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
            -DBUILD_SHARED_LIBS=ON
            -DENABLE_NETCDF_4=ON
            -DENABLE_DAP=ON
            -DENABLE_TESTS=OFF
            -DBUILD_UTILITIES=OFF
            -DHDF5_ROOT=${HDF5_INSTALL_DIR}
            -DHDF5_DIR=${HDF5_INSTALL_DIR}/lib/cmake/hdf5
            -DCMAKE_PREFIX_PATH=${HDF5_INSTALL_DIR}
            # Suppress compiler warnings (optional)
            -DCMAKE_C_FLAGS=-w
        BUILD_COMMAND ${CMAKE_COMMAND} --build . --parallel
        INSTALL_COMMAND ${CMAKE_COMMAND} --install .
    )

    # Create imported targets for easier linking
    add_library(HDF5::HDF5 SHARED IMPORTED)
    set_target_properties(HDF5::HDF5 PROPERTIES
        IMPORTED_LOCATION ${HDF5_INSTALL_DIR}/lib/libhdf5${CMAKE_SHARED_LIBRARY_SUFFIX}
        INTERFACE_INCLUDE_DIRECTORIES ${HDF5_INSTALL_DIR}/include
    )

    add_library(NetCDF::NetCDF SHARED IMPORTED)
    set_target_properties(NetCDF::NetCDF PROPERTIES
        IMPORTED_LOCATION ${NETCDF_C_INSTALL_DIR}/lib/libnetcdf${CMAKE_SHARED_LIBRARY_SUFFIX}
        INTERFACE_INCLUDE_DIRECTORIES ${NETCDF_C_INSTALL_DIR}/include
        INTERFACE_LINK_LIBRARIES HDF5::HDF5
    )

    # Add dependencies to ensure external projects are built before imported targets are used
    add_dependencies(HDF5::HDF5 HDF5)
    add_dependencies(NetCDF::NetCDF NetCDFC)

    # Set variables for compatibility (if your code expects these)
    set(HDF5_INCLUDE_DIR ${HDF5_INSTALL_DIR}/include)
    set(HDF5_LIBRARY ${HDF5_INSTALL_DIR}/lib/libhdf5${CMAKE_SHARED_LIBRARY_SUFFIX})
    set(NETCDF_C_INCLUDE_DIR ${NETCDF_C_INSTALL_DIR}/include)
    set(NETCDF_C_LIBRARY ${NETCDF_C_INSTALL_DIR}/lib/libnetcdf${CMAKE_SHARED_LIBRARY_SUFFIX})

    # Add compile definitions
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSTORMM_USE_NETCDF")
    if(STORMM_ENABLE_CUDA)
        set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -DSTORMM_USE_NETCDF")
    endif()
endif()

# --- RDKit Integration ---
if (STORMM_ENABLE_RDKIT)
    message("-- RDKit support will be enabled")

    find_package(RDKit REQUIRED)

    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -DSTORMM_USE_RDKIT -DSTORMM_USE_RDKIT")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSTORMM_USE_RDKIT -DSTORMM_USE_RDKIT")
endif()

#--- Define any global settings / flags ---
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)  # set the -fPIC flag

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -O0 -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3")
if (CMAKE_SHOW_FALLTHROUGH)
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -Wimplicit-fallthrough")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -Wimplicit-fallthrough")
endif()
set(STORMM_VERSION "0.4")

#--- Functions for detecting a program's availability and version, e.g. gcc version ---
set(PROGRAM_NAME my_program)

function(${PROGRAM_NAME}_version output_var)

  # Run the command to get the pluginval version
  execute_process(COMMAND ${PROGRAM_NAME} --version
                  OUTPUT_VARIABLE ${PROGRAM_NAME}_VERSION_RAW_OUTPUT)

  # Extracts the version from the output of the command run before
  string(STRIP ${${PROGRAM_NAME}_VERSION_OUTPUT} ${PROGRAM_NAME}_VERSION_OUTPUT)

  set(${output_var} ${${PROGRAM_NAME}_VERSION_OUTPUT} PARENT_SCOPE)
endfunction()

# Create the memory check script
file(WRITE ${CMAKE_BINARY_DIR}/check_memory.cmake
"if(CMAKE_HOST_SYSTEM_NAME STREQUAL \"Linux\")
    execute_process(
        COMMAND grep MemTotal /proc/meminfo
        OUTPUT_VARIABLE TOTAL_MEM_KB
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    execute_process(
        COMMAND grep MemAvailable /proc/meminfo
        OUTPUT_VARIABLE AVAILABLE_MEM_KB
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    string(REGEX REPLACE \"[^0-9]\" \"\" TOTAL_MEM_KB \"\${TOTAL_MEM_KB}\")
    string(REGEX REPLACE \"[^0-9]\" \"\" AVAILABLE_MEM_KB \"\${AVAILABLE_MEM_KB}\")

    math(EXPR USED_PERCENT \"100 - (\${AVAILABLE_MEM_KB} * 100 / \${TOTAL_MEM_KB})\")

    if(USED_PERCENT GREATER 70)
        message(WARNING \"System memory usage is at \${USED_PERCENT}% (above 70% threshold)\")
        message(STATUS \"Available memory: \${AVAILABLE_MEM_KB}KB of \${TOTAL_MEM_KB}KB total\")
        message(FATAL_ERROR \"Memory usage too high. Please free up memory or use -DSTORMM_IGNORE_RAM_LIMITS=ON to override\")
    else()
        message(STATUS \"Memory usage at \${USED_PERCENT}% - proceeding with compilation\")
    endif()
elseif(CMAKE_HOST_SYSTEM_NAME STREQUAL \"Darwin\")
    execute_process(
        COMMAND sh -c \"vm_stat | grep 'Pages free:' | awk '{print \$3}'\"
        OUTPUT_VARIABLE FREE_PAGES
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    execute_process(
        COMMAND sh -c \"vm_stat | grep 'Pages inactive:' | awk '{print \$3}'\"
        OUTPUT_VARIABLE INACTIVE_PAGES
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    execute_process(
        COMMAND sysctl -n hw.memsize
        OUTPUT_VARIABLE TOTAL_MEM_B
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    string(REGEX REPLACE \"[^0-9]\" \"\" FREE_PAGES \"\${FREE_PAGES}\")
    string(REGEX REPLACE \"[^0-9]\" \"\" INACTIVE_PAGES \"\${INACTIVE_PAGES}\")

    math(EXPR AVAILABLE_MEM_B \"(\${FREE_PAGES} + \${INACTIVE_PAGES}) * 4096\")
    math(EXPR AVAILABLE_PERCENT \"(\${AVAILABLE_MEM_B} * 100 / \${TOTAL_MEM_B})\")

    if(AVAILABLE_PERCENT GREATER 70)
        message(WARNING \"System memory availability is at \${AVAILABLE_PERCENT}% (above 70% threshold)\")
        message(STATUS \"Available memory: \${AVAILABLE_MEM_B} bytes of \${TOTAL_MEM_B} bytes total\")
        message(FATAL_ERROR \"Memory availability too low. Please free up memory or use -DSTORMM_IGNORE_RAM_LIMITS=ON to override\")
    else()
        message(STATUS \"Memory availability at \${AVAILABLE_PERCENT}% - proceeding with compilation\")
        message(STATUS \"Available memory: \${AVAILABLE_MEM_B} bytes of \${TOTAL_MEM_B} bytes total\")
    endif()
endif()
")

if (STORMM_ENABLE_CUDA)
    message("-- CUDA support will be enabled")

    enable_language(CUDA)

    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)

    # Enable split compilation to reduce memory usage
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -DSTORMM_USE_HPC -DSTORMM_USE_CUDA")
    set(CMAKE_CUDA_FLAGS_DEBUG "${CMAKE_CUDA_FLAGS_DEBUG} --compiler-options=-O0 -g")

    # Detect GPU architecture using nvidia-smi
    execute_process(
        COMMAND nvidia-smi --query-gpu=compute_cap --format=csv,noheader,nounits
        OUTPUT_VARIABLE GPU_COMPUTE_CAP_RAW
        OUTPUT_STRIP_TRAILING_WHITESPACE
        RESULT_VARIABLE NVIDIA_SMI_RESULT
    )
    if (NOT DEFINED CUSTOM_GPU_ARCH)
      if (NVIDIA_SMI_RESULT EQUAL 0 AND GPU_COMPUTE_CAP_RAW)
            # Split the output into a list based on newlines
            string(REPLACE "\n" ";" GPU_COMPUTE_CAP_LIST "${GPU_COMPUTE_CAP_RAW}")
            # Create a list to store processed architectures
            set(ARCH_LIST "")

            # Process each compute capability and remove .
            foreach(CAP IN LISTS GPU_COMPUTE_CAP_LIST)
                  string(REPLACE "." "" GPU_ARCH "${CAP}")
                  list(APPEND ARCH_LIST ${GPU_ARCH})
            endforeach()

            list(REMOVE_DUPLICATES ARCH_LIST)

            # Set the architectures in the CMake variable
            set(CMAKE_CUDA_ARCHITECTURES ${ARCH_LIST})
            message("-- Detected GPU Architecture(s): ${GPU_COMPUTE_CAP_RAW} (compute_${ARCH_LIST})")
        else()
            # Fallback to default if detection fails
            set(CMAKE_CUDA_ARCHITECTURES "52")
            message("-- Could not detect GPU architecture, defaulting to compute capability 5.2")
            message("   Set CUSTOM_GPU_ARCH if you need a different architecture")
        endif()
    else()
        set(CMAKE_CUDA_ARCHITECTURES "${CUSTOM_GPU_ARCH}")
        message("-- Using user-specified GPU architecture: compute_${CUSTOM_GPU_ARCH}")
    endif()

    if (NOT DEFINED CUSTOM_NVCC_THREADS)
        set(CMAKE_NVCC_THREADS 1)
    else()
        set(CMAKE_NVCC_THREADS ${CUSTOM_NVCC_THREADS})
    endif()

    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --threads ${CMAKE_NVCC_THREADS}")
    set(CMAKE_CUDA_FLAGS_DEBUG "${CMAKE_CUDA_FLAGS_DEBUG} --threads ${CMAKE_NVCC_THREADS}")
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --use_fast_math")

    if (CMAKE_SHOW_FALLTHROUGH)
        set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Wimplicit-fallthrough")
    endif()
    if (CMAKE_SHOW_PTXAS)
        set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --ptxas-options=\"-v\"")
    endif()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSTORMM_USE_HPC -DSTORMM_USE_CUDA")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSTORMM_USE_HPC -DSTORMM_USE_CUDA")

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
find_package(Git)

if (GIT_EXECUTABLE)
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

if (NOT DEFINED STORMM_VERSION)
    set(STORMM_VERSION 0+unknown)
    message(WARNING "Could not find STORMM_VERSION using git... using default version \"${STORMM_VERSION}\".")
endif()

configure_file(
        "${CMAKE_SOURCE_DIR}/src/version.h.in"
        "${CMAKE_BINARY_DIR}/version.h"
        @ONLY)

# --- Source Collection and Build Phases ---
# First define all sources, separating CPU and GPU code for phased compilation.
# This separation enables memory optimization by allowing CPU compilation to
# complete and free memory before GPU compilation begins.

# Source file collection
# Using CONFIGURE_DEPENDS tells CMake to check for new/removed source files
file(GLOB_RECURSE STORMM_CPU_SOURCES
    CONFIGURE_DEPENDS
    "${CMAKE_SOURCE_DIR}/src/*.c"    # C source files
    "${CMAKE_SOURCE_DIR}/src/*.cpp"  # C++ source files
    "${CMAKE_SOURCE_DIR}/src/*.tpp"  # Template implementations
)

# Collect GPU sources if CUDA is enabled
if(STORMM_ENABLE_CUDA)
    file(GLOB_RECURSE STORMM_CUDA_SOURCES
        CONFIGURE_DEPENDS
        "${CMAKE_SOURCE_DIR}/src/*.cu"   # CUDA source files
        "${CMAKE_SOURCE_DIR}/src/*.cuh"  # CUDA headers
        "${CMAKE_SOURCE_DIR}/src/*.cui"  # CUDA inline implementations
    )
endif()

# --- CPU Compilation Phase ---
# Create object library for CPU sources to compile them independently.
# Using object libraries allows us to:
# 1. Compile CPU code first and completely
# 2. Free memory after CPU compilation
# 3. Reuse the compiled objects in the final library
add_library(${PROJECT_NAME}_cpu_objects OBJECT ${STORMM_CPU_SOURCES})
target_include_directories(${PROJECT_NAME}_cpu_objects
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:include/${PROJECT_NAME}>
)

if(STORMM_ENABLE_RDKIT)
    target_link_libraries(${PROJECT_NAME}_cpu_objects PUBLIC RDKit::SmilesParse)
endif()

# Mark CPU compilation completion
# This target serves as a checkpoint between CPU and GPU compilation phases
add_custom_target(${PROJECT_NAME}_cpu_complete
    DEPENDS ${PROJECT_NAME}_cpu_objects
    COMMAND ${CMAKE_COMMAND} -E echo "CPU compilation complete. Memory can be freed."
)

# Modify the memory check targets to only run if RAM limits aren't ignored
if(NOT STORMM_IGNORE_RAM_LIMITS)
    add_custom_target(${PROJECT_NAME}_check_memory_cpu
        COMMAND ${CMAKE_COMMAND} -P ${CMAKE_BINARY_DIR}/check_memory.cmake
    )
    add_dependencies(${PROJECT_NAME}_cpu_objects ${PROJECT_NAME}_check_memory_cpu)

    if(STORMM_ENABLE_CUDA)
        add_custom_target(${PROJECT_NAME}_check_memory_gpu
            COMMAND ${CMAKE_COMMAND} -P ${CMAKE_BINARY_DIR}/check_memory.cmake
        )
        add_dependencies(${PROJECT_NAME}_cuda_objects ${PROJECT_NAME}_check_memory_gpu)
    endif()
endif()

# --- Main Library Target ---
# This will contain both CPU and GPU components when fully linked
add_library(${PROJECT_NAME} SHARED)

# --- GPU Compilation Phase ---
if(STORMM_ENABLE_CUDA)
    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)

    # Create CUDA object library
    # This compilation only starts after CPU phase is complete
    add_library(${PROJECT_NAME}_cuda_objects OBJECT ${STORMM_CUDA_SOURCES})
    target_include_directories(${PROJECT_NAME}_cuda_objects
        PUBLIC
        ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
    )

    # Ensure GPU compilation waits for CPU completion
    # This dependency makes 'make -j' compile CPU first, then GPU
    add_dependencies(${PROJECT_NAME}_cuda_objects ${PROJECT_NAME}_cpu_complete)

    # Mark GPU compilation completion
    add_custom_target(${PROJECT_NAME}_cuda_complete
        DEPENDS ${PROJECT_NAME}_cuda_objects
        COMMAND ${CMAKE_COMMAND} -E echo "CUDA compilation complete. Memory can be freed."
    )

    # Combine all objects into final library
    target_sources(${PROJECT_NAME}
        PRIVATE
        $<TARGET_OBJECTS:${PROJECT_NAME}_cpu_objects>
        $<TARGET_OBJECTS:${PROJECT_NAME}_cuda_objects>
    )

    # Ensure final linking waits for all compilation
    add_dependencies(${PROJECT_NAME} ${PROJECT_NAME}_cuda_complete)

    target_link_libraries(${PROJECT_NAME} PUBLIC ${STORMM_CUDA_LIBS})
else()
    # CPU-only build
    target_sources(${PROJECT_NAME}
        PRIVATE
        $<TARGET_OBJECTS:${PROJECT_NAME}_cpu_objects>
    )
    add_dependencies(${PROJECT_NAME} ${PROJECT_NAME}_cpu_complete)
endif()

# Add include directories to main library
target_include_directories(${PROJECT_NAME}
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:include/${PROJECT_NAME}>
)

if(STORMM_ENABLE_CUDA)
    target_include_directories(${PROJECT_NAME}
        PUBLIC
        ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}
    )
endif()

# Installation rules
install(TARGETS ${PROJECT_NAME}
    EXPORT ${PROJECT_NAME_UPPER}Targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    BUNDLE DESTINATION ${CMAKE_INSTALL_BINDIR}
)

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

#--- Build the Python bridge ---
if (STORMM_ENABLE_PYTHON)
    # Build a small C ABI layer for Python ctypes clients.  This target links against the main
    # STORMM shared library and exports stable, pointer-based entry points.
    add_library(stormm_pybridge SHARED python/cpp/stormm_pybridge.cpp)
    target_include_directories(stormm_pybridge
        PRIVATE
        ${CMAKE_SOURCE_DIR}/src
    )
    target_link_libraries(stormm_pybridge
        PRIVATE
        ${PROJECT_NAME}
    )
    set_target_properties(stormm_pybridge PROPERTIES
        OUTPUT_NAME "stormm_pybridge"
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
    )
endif()
