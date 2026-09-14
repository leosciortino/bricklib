VSTARGET(N3AXIS ../stencils/3axis.cpp ../stencils/3axis-out.cpp)
add_executable(single-cpu cpu.cpp ${VSTARGET_N3AXIS_OUTPUT} ../stencils/copy.cpp)
target_link_libraries(single-cpu brickhelper)
set_target_properties(single-cpu
        PROPERTIES OUTPUT_NAME cpu)

add_executable(single-mpi mpi.cpp ../src/brick-mpi.cpp)
target_link_libraries(single-mpi brick brickhelper ${MPI_CXX_LIBRARIES})
set_target_properties(single-mpi
        PROPERTIES OUTPUT_NAME mpi)

if (HAS_SYCL)
    set(SYCL3AXIS_INCLUDE_DIRS
            "${SYCL_DIR}/include"
            "${SYCL_DIR}/include/sycl")
    set(SYCL3AXIS_COMPILE_OPTIONS
            "-fsycl")
    VSTARGET(SYCL3AXIS ../stencils/3axis_sycl.cpp ../stencils/3axis_sycl-out.cpp)
    add_executable(single-sycl sycl.cpp ${VSTARGET_SYCL3AXIS_OUTPUT})
    target_link_libraries(single-sycl brick brickhelper stdc++ sycl)
    target_include_directories(single-sycl 
            PUBLIC "${SYCL_DIR}/include"
            PUBLIC "${SYCL_DIR}/include/sycl"
            PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/stencils")
    target_link_directories(single-sycl PUBLIC "${SYCL_DIR}/lib")
    # sycl need option for both phase
    target_compile_options(single-sycl PUBLIC "-fsycl")
    target_link_options(single-sycl PUBLIC "-fsycl")
    set_target_properties(single-sycl
            PROPERTIES OUTPUT_NAME sycl)
endif ()

if (OpenCL_FOUND)
    add_executable(single-opencl opencl.cpp ../stencils/3axis_cl.cpp)
    target_compile_definitions(single-opencl PRIVATE "CL_KRNL_OPTIONS=\"-I${PROJECT_SOURCE_DIR}/include\"")
    VSTARGET(CL3AXIS ../stencils/3axis_cl_krnl.c ${CMAKE_CURRENT_BINARY_DIR}/3axis_cl_krnl-out.c)
    add_custom_target(cl3axis
            SOURCES ${VSTARGET_CL3AXIS_OUTPUT})
    target_link_libraries(single-opencl brick brickhelper ${OpenCL_LIBRARIES})
    target_include_directories(single-opencl PUBLIC "${OpenCL_INCLUDE_DIRS}")
    add_dependencies(single-opencl cl3axis)
    set_target_properties(single-opencl
            PROPERTIES OUTPUT_NAME opencl)
endif ()

if (CMAKE_CUDA_COMPILER)
    set(C3AXIS_INCLUDE_DIRS
            "${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}"
            )
    VSTARGET(C3AXIS ../stencils/3axis.cu ../stencils/3axis-out.cu)
    add_executable(single-cuda cuda.cpp ${VSTARGET_C3AXIS_OUTPUT})
    target_link_libraries(single-cuda brickhelper)
    set_target_properties(single-cuda
            PROPERTIES OUTPUT_NAME cuda)
    target_include_directories(single-cuda PUBLIC
            "${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}"
            )

        # New sparse CUDA benchmark
    set(SPARSE_C3AXIS_INCLUDE_DIRS
        "${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}"
    )
    set(SPARSE_C3AXIS_DEPENDS
        "${CMAKE_CURRENT_SOURCE_DIR}/../stencils/python/star_7.py"
        "${CMAKE_CURRENT_SOURCE_DIR}/../stencils/python/star_13.py"
        "${CMAKE_CURRENT_SOURCE_DIR}/../stencils/python/star_19.py"
        "${CMAKE_CURRENT_SOURCE_DIR}/../stencils/python/star_25.py"
        "${CMAKE_CURRENT_SOURCE_DIR}/../stencils/python/cube_27.py"
        "${CMAKE_CURRENT_SOURCE_DIR}/../stencils/python/cube_125.py"
    )

    VSTARGET(
        SPARSE_C3AXIS
        ../stencils/3axis_sparse.cu
        ../stencils/3axis_sparse-out.cu
    )

    add_executable(
        single-sparse-cuda
        sparse_cuda.cpp
        ${VSTARGET_SPARSE_C3AXIS_OUTPUT}
    )
    set_source_files_properties(sparse_cuda.cpp PROPERTIES LANGUAGE CUDA)

    target_link_libraries(
        single-sparse-cuda
        brickhelper
    )
    target_link_options(single-sparse-cuda PRIVATE
        "-L${CMAKE_CUDA_COMPILER_TOOLKIT_ROOT}/lib64"
    )

    set_target_properties(
        single-sparse-cuda
        PROPERTIES OUTPUT_NAME sparse-cuda
    )

    target_include_directories(
        single-sparse-cuda PUBLIC
        "${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}"
    )
endif ()

if (BRICK_USE_HIP)
    set(H3AXIS_INCLUDE_DIRS "${HIP_PATH}/include")
    VSTARGET(H3AXIS ../stencils/3axis.hip.cpp ../stencils/3axis-out.hip.cpp)
    set(HIP_SRCS hip.cpp ${VSTARGET_H3AXIS_OUTPUT})
    set_source_files_properties(${HIP_SRCS} PROPERTIES HIP_SOURCE_PROPERTY_FORMAT 1)
    hip_add_executable(single-hip ${HIP_SRCS} CLANG_OPTIONS ${CMAKE_CXX_FLAGS})
    target_link_libraries(single-hip brickhelper hip::host)
    set_target_properties(single-hip
            PROPERTIES OUTPUT_NAME hip)

    set(SPARSE_H3AXIS_INCLUDE_DIRS "${HIP_PATH}/include")
    VSTARGET(SPARSE_H3AXIS
            ../stencils/3axis_sparse.hip.cpp
            ../stencils/3axis_sparse-out.hip.cpp)
    set(SPARSE_HIP_SRCS sparse_hip.cpp ${VSTARGET_SPARSE_H3AXIS_OUTPUT})
    set_source_files_properties(${SPARSE_HIP_SRCS}
            PROPERTIES HIP_SOURCE_PROPERTY_FORMAT 1)
    hip_add_executable(single-sparse-hip ${SPARSE_HIP_SRCS}
            HIPCC_OPTIONS --offload-arch=gfx942
            CLANG_OPTIONS ${CMAKE_CXX_FLAGS})
    target_link_libraries(single-sparse-hip brickhelper hip::host)
    set_target_properties(single-sparse-hip
            PROPERTIES OUTPUT_NAME sparse-hip)
endif ()

if (USE_DPC)
   VSTARGET(D3AXIS ../stencils/3axis_dpc.cpp ../stencils/3axis_dpc-out.cpp)
   set(CMAKE_CXX_COMPILER "dpcpp")
   set(CMAKE_LINKER "dpcpp") 
   set(CMAKE_AR "llvm-ar") 
   set(CMAKE_RANLIB "llvm-ranlib") 
   set(CMAKE_CXX_STANDARD 17) 
   add_executable(single-dpc dpc.cpp ${VSTARGET_D3AXIS_OUTPUT})
   target_link_libraries(single-dpc brickhelper)
   target_include_directories(single-dpc PUBLIC "$(which dpcpp)/../../include/sycl/")
   set_target_properties(single-dpc PROPERTIES OUTPUT_NAME dpc)
endif ()
