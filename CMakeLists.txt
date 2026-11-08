cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
# 需要 C 语言：extern/catima 与 extern/fftw 都是 C/C++ 源码
project(SpyralPP VERSION 1.0 LANGUAGES C CXX)

#----------------------------------------------------------------------
# 基本设置
#----------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# 安装后程序用 $ORIGIN/../lib 找 librk4.so
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib")

include(GNUInstallDirs)

#----------------------------------------------------------------------
# 依赖
#   唯一需要预装的是 ROOT（含 Minuit）
#   其余全部自带源码，位于 extern/：
#     extern/catima       catima 1.8.0          能量损失 / 射程
#     extern/fftw         FFTW 3.3.10           波形滤波（双精度）
#     extern/hdf5         HDF5 1.12.0           预编译静态库 + 头文件
#     extern/highfive     HighFive              HDF5 的 C++ 封装（header-only）
#     extern/rk4_solver   RK4 积分器             header-only
#----------------------------------------------------------------------
find_package(ROOT REQUIRED)
find_package(Threads REQUIRED)

find_library(ROOT_MINUIT_LIB NAMES Minuit HINTS ${ROOT_LIBRARY_DIR})
find_library(ROOT_MINUIT2_LIB NAMES Minuit2 HINTS ${ROOT_LIBRARY_DIR})
if(NOT ROOT_MINUIT2_LIB)
    message(FATAL_ERROR "ROOT Minuit2 is required for independent parallel trajectory fits")
endif()
if(NOT ROOT_MINUIT_LIB)
    message(FATAL_ERROR "找不到 libMinuit（ROOT 的 TMinuit 不在默认链接列表里）")
endif()

#----------------------------------------------------------------------
# 内置第三方库（源码形式）
#   两个子项目都会自己 declare project() 并带 install 规则，
#   这里统一按静态库编译，最后一起塞进 librk4.so。
#----------------------------------------------------------------------
set(BUILD_SHARED_LIBS OFF CACHE BOOL "内置第三方库一律编译成静态库" FORCE)

# --- catima ---
set(_catima_dir "${CMAKE_CURRENT_SOURCE_DIR}/extern/catima")
if(NOT EXISTS "${_catima_dir}/CMakeLists.txt")
    message(FATAL_ERROR "extern/catima 源码不存在")
endif()
set(TESTS        OFF CACHE BOOL "" FORCE)
set(EXAMPLES     OFF CACHE BOOL "" FORCE)
set(APPS         OFF CACHE BOOL "" FORCE)
set(GLOBAL       OFF CACHE BOOL "" FORCE)
set(GENERATE_DATA OFF CACHE BOOL "" FORCE)
add_subdirectory("${_catima_dir}" "${CMAKE_BINARY_DIR}/extern/catima" EXCLUDE_FROM_ALL)

# --- FFTW ---
set(_fftw_dir "${CMAKE_CURRENT_SOURCE_DIR}/extern/fftw")
if(NOT EXISTS "${_fftw_dir}/CMakeLists.txt")
    message(FATAL_ERROR "extern/fftw 源码不存在")
endif()
set(BUILD_TESTS        OFF CACHE BOOL "" FORCE)
set(ENABLE_FLOAT       OFF CACHE BOOL "" FORCE)
set(ENABLE_LONG_DOUBLE OFF CACHE BOOL "" FORCE)
set(ENABLE_QUAD_PRECISION OFF CACHE BOOL "" FORCE)
set(ENABLE_THREADS     OFF CACHE BOOL "" FORCE)
set(ENABLE_OPENMP      OFF CACHE BOOL "" FORCE)
set(WITH_COMBINED_THREADS OFF CACHE BOOL "" FORCE)
set(DISABLE_FORTRAN    ON  CACHE BOOL "" FORCE)
set(ENABLE_SSE2        ON  CACHE BOOL "" FORCE)
set(ENABLE_AVX         ON  CACHE BOOL "" FORCE)
set(ENABLE_AVX2        ON  CACHE BOOL "" FORCE)
add_subdirectory("${_fftw_dir}" "${CMAKE_BINARY_DIR}/extern/fftw" EXCLUDE_FROM_ALL)

# FFTW 的 CMakeLists 只给了 INSTALL_INTERFACE 的头文件路径，
# 源码树内链接的调用方需要手动补一个 api/ 目录。
target_include_directories(fftw3 INTERFACE $<BUILD_INTERFACE:${_fftw_dir}/api>)

message(STATUS "ROOT      : ${ROOT_DIR}")
message(STATUS "Minuit    : ${ROOT_MINUIT_LIB} ${ROOT_MINUIT2_LIB}")
message(STATUS "catima    : 内置源码 extern/catima")
message(STATUS "FFTW      : 内置源码 extern/fftw")

#----------------------------------------------------------------------
# 选项
#----------------------------------------------------------------------
option(SPYRALPP_WITH_HDF5 "启用 HDF5 读写（内置 extern/hdf5）" ON)
option(SPYRALPP_WITH_LEGACY_ANALYSIS "编译 ParticleDynamics 的遗留分析代码" OFF)
option(SPYRALPP_BUILD_APPS "构建命令行程序" ON)
option(SPYRALPP_BUILD_EXPERIMENTAL "构建研究用实验程序（需要 experimental/ 源码）" OFF)

#----------------------------------------------------------------------
# 核心库
#----------------------------------------------------------------------
add_library(rk4 SHARED
    src/core/MassData.C
    src/core/DataPath.C
    src/core/Reaction.C
    src/io/TreeIO.C
    src/io/CloudIO.C
    src/io/EventMetaIO.C
    src/io/WorkspaceIO.C
    src/io/WorkspaceH5.C
    src/io/Config.C
    src/io/ResultIO.C
    src/io/StageIO.C
    src/sim/TargetDensity.C
    src/sim/Stopping.C
    src/sim/KinematicsGenerator.C
    src/sim/PadPlane.C
    src/sim/Response.C
    src/sim/Simulator.C
    src/sim/ParticleDynamics.C
    src/reco/WaveProcessor.C
    src/reco/PulseShape.C
    src/reco/PeakFinder.C
    src/reco/TraceReader.C
    src/reco/H5TraceReader.C
    src/reco/PointCloudBuilder.C
    src/reco/Clusterer.C
    src/reco/ClustererSpyral.C
    src/reco/Hdbscan.C
    src/reco/Geometry.C
    src/reco/Spline.C
    src/reco/EstimatorReinsch.C
    src/reco/Estimator.C
    src/reco/ParticleId.C
    src/reco/TrackMesh.C
    src/reco/TrackFitter.C
    src/reco/KalmanFitter.C
)

target_include_directories(rk4 PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${ROOT_INCLUDE_DIRS}
    # Kalman 滤波：OpenKF 的 UKF 工具 + 它依赖的 Eigen（都是源码）
    ${CMAKE_CURRENT_SOURCE_DIR}/extern/openkf
    ${CMAKE_CURRENT_SOURCE_DIR}/extern/eigen
)

# rk4_solver's default Integrator stores each fixed-size RK4 work array using
# raw `new[]` references without a destructor. UKF path prediction constructs
# many short-lived integrators, so the default leaks heap on every RK step.
# The state dimension is six here; keep its small work arrays on the stack.
target_compile_definitions(rk4 PUBLIC DO_NOT_USE_HEAP)

target_link_libraries(rk4 PUBLIC
    ${ROOT_LIBRARIES}
    ${ROOT_MINUIT_LIB} ${ROOT_MINUIT2_LIB}
    catima::catima
    fftw3
    Threads::Threads
    m
)

if(SPYRALPP_WITH_LEGACY_ANALYSIS)
    target_compile_definitions(rk4 PUBLIC RK4_WITH_LEGACY_ANALYSIS)
    message(STATUS "Legacy    : 开启")
else()
    message(STATUS "Legacy    : 关闭")
endif()

#----------------------------------------------------------------------
# HDF5（内置源码）
#----------------------------------------------------------------------
if(SPYRALPP_WITH_HDF5)
    set(_h5_inc "${CMAKE_CURRENT_SOURCE_DIR}/extern/hdf5/include")
    set(_h5_lib "${CMAKE_CURRENT_SOURCE_DIR}/extern/hdf5/lib/libhdf5.a")
    if(NOT EXISTS "${_h5_inc}/hdf5.h" OR NOT EXISTS "${_h5_lib}")
        find_package(HDF5 REQUIRED COMPONENTS C)
        set(_h5_inc ${HDF5_INCLUDE_DIRS})
        set(_h5_lib ${HDF5_LIBRARIES})
    endif()

    target_sources(rk4 PRIVATE src/io/H5Writer.C)
    target_include_directories(rk4 PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/extern/highfive/include
        ${_h5_inc})
    target_link_libraries(rk4 PUBLIC ${_h5_lib})
    target_compile_definitions(rk4 PUBLIC RK4_WITH_HDF5)
    message(STATUS "HDF5      : 开启（${_h5_lib}）")
else()
    target_sources(rk4 PRIVATE src/io/H5WriterStub.C)
    message(STATUS "HDF5      : 关闭")
endif()

#----------------------------------------------------------------------
# 命令行程序
#----------------------------------------------------------------------
if(SPYRALPP_BUILD_APPS)
    foreach(app simulate analyze waveform csvtrace merge estimate_csv pid workspace)
        add_executable(rk4_${app} src/app/${app}.cpp)
        target_link_libraries(rk4_${app} PRIVATE rk4)
        install(TARGETS rk4_${app} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
    endforeach()
    if(SPYRALPP_BUILD_EXPERIMENTAL)
    add_executable(experimental_minuit_objective_order_probe
        experimental/minuit_multiturn_audit/objective_order_probe.cpp)
    target_include_directories(experimental_minuit_objective_order_probe PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
    target_link_libraries(experimental_minuit_objective_order_probe PRIVATE rk4)
    add_executable(experimental_minuit_objective_order_batch_probe
        experimental/minuit_multiturn_audit/objective_order_batch_probe.cpp)
    target_include_directories(experimental_minuit_objective_order_batch_probe PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
    target_link_libraries(experimental_minuit_objective_order_batch_probe PRIVATE rk4)
    add_executable(experimental_minuit_direct_rk4
        experimental/minuit_multiturn_audit/direct_minuit/direct_rk4_minuit.cpp)
    target_include_directories(experimental_minuit_direct_rk4 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
    target_link_libraries(experimental_minuit_direct_rk4 PRIVATE rk4)
    add_executable(experimental_ukf_cluster_snapshot
        experimental/ukf_cluster_merging/cluster_stage_adapter/export_snapshot.cpp)
    target_include_directories(experimental_ukf_cluster_snapshot PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
    target_link_libraries(experimental_ukf_cluster_snapshot PRIVATE rk4)
    add_executable(experimental_phase1a_replay_finalize
        experimental/cluster_recovery/phase1a/replay_finalize.cpp)
    target_include_directories(experimental_phase1a_replay_finalize PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
    target_link_libraries(experimental_phase1a_replay_finalize PRIVATE rk4)
    add_executable(zpmt_ukf_edge_validator
        experimental/zsliced_multitarget/ukf_edge_validator.cpp)
    target_include_directories(zpmt_ukf_edge_validator PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
    target_link_libraries(zpmt_ukf_edge_validator PRIVATE rk4)
    add_executable(phase1c_c1p_validator
        experimental/multi_model_track_reconstruction/c1p_validator.cpp)
    target_include_directories(phase1c_c1p_validator PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
    target_link_libraries(phase1c_c1p_validator PRIVATE rk4)
    add_executable(experimental_ukf_fragmentation
        experimental/ukf_cluster_merging/artificial_fragmentation/run_original_baselines.cpp)
    target_include_directories(experimental_ukf_fragmentation PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
    target_link_libraries(experimental_ukf_fragmentation PRIVATE rk4)
    add_executable(experimental_ukf_gap_deletion
        experimental/ukf_cluster_merging/artificial_fragmentation/run_gap_deletion.cpp)
    target_include_directories(experimental_ukf_gap_deletion PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
    target_link_libraries(experimental_ukf_gap_deletion PRIVATE rk4)
    add_executable(experimental_ukf_merge_case
        experimental/ukf_cluster_merging/ukf_merger/run_ukf_case.cpp
        experimental/ukf_cluster_merging/ukf_merger/UkfClusterMerger.cpp)
    target_include_directories(experimental_ukf_merge_case PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/experimental/ukf_cluster_merging/ukf_merger)
    target_link_libraries(experimental_ukf_merge_case PRIVATE rk4)
    add_executable(experimental_ukf_full_batch
        experimental/ukf_cluster_merging/ukf_merger/run_ukf_full_batch.cpp
        experimental/ukf_cluster_merging/ukf_merger/UkfClusterMerger.cpp)
    target_include_directories(experimental_ukf_full_batch PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/experimental/ukf_cluster_merging/ukf_merger)
    target_link_libraries(experimental_ukf_full_batch PRIVATE rk4)
    add_executable(experimental_merge_edge_replay
        experimental/cluster_merge_truth_audit/merge_edge_ablation/replay_one_edge.cpp)
    target_include_directories(experimental_merge_edge_replay PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
    target_link_libraries(experimental_merge_edge_replay PRIVATE rk4)
    endif()
endif()

#----------------------------------------------------------------------
# 安装
#----------------------------------------------------------------------
include(GNUInstallDirs)

install(TARGETS rk4
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(DIRECTORY src/core src/sim src/reco src/io
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/spyralpp
        FILES_MATCHING PATTERN "*.h")

install(DIRECTORY data
        DESTINATION ${CMAKE_INSTALL_DATADIR}/spyralpp
        FILES_MATCHING PATTERN "*.csv" PATTERN "*.txt" PATTERN "amdc_2020")

install(DIRECTORY configs
        DESTINATION ${CMAKE_INSTALL_DATADIR}/spyralpp
        FILES_MATCHING PATTERN "*.cfg" PATTERN "*.txt")

# 一键脚本（跑批量任务用）：<prefix>/share/spyralpp/scripts/
install(DIRECTORY scripts
        DESTINATION ${CMAKE_INSTALL_DATADIR}/spyralpp
        FILES_MATCHING PATTERN "*.sh"
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                    GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

include(CTest)
if(BUILD_TESTING)
    add_executable(test_mesh_context test/testMeshContext.cpp)
    target_link_libraries(test_mesh_context PRIVATE rk4)
    add_test(NAME mesh_context COMMAND test_mesh_context)
    set_tests_properties(mesh_context PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif()
if(BUILD_TESTING)
    add_executable(test_threaded_stopping test/testThreadedStopping.cpp)
    target_link_libraries(test_threaded_stopping PRIVATE rk4)
    add_test(NAME threaded_stopping COMMAND test_threaded_stopping)
    set_tests_properties(threaded_stopping PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif()
if(BUILD_TESTING)
    foreach(check cloud_validation lifecycle)
        if(check STREQUAL "cloud_validation")
            set(check_source test/testCloudValidation.cpp)
        else()
            set(check_source test/testLifecycle.cpp)
        endif()
        add_executable(test_${check} ${check_source})
        target_link_libraries(test_${check} PRIVATE rk4)
        add_test(NAME ${check} COMMAND test_${check})
        set_tests_properties(${check} PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    endforeach()
endif()
if(BUILD_TESTING)
    add_executable(test_target_mixtures test/testTargetMixtures.cpp)
    target_link_libraries(test_target_mixtures PRIVATE rk4)
    add_test(NAME target_mixtures COMMAND test_target_mixtures)
    set_tests_properties(target_mixtures PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif()
if(BUILD_TESTING)
    add_executable(test_catima_spline test/testCatimaSpline.cpp)
    target_link_libraries(test_catima_spline PRIVATE rk4)
    add_test(NAME catima_spline COMMAND test_catima_spline)
endif()
if(BUILD_TESTING)
    add_executable(test_kinematics_regression test/kinematics_regression.cpp)
    add_executable(test_detector_regression test/detector_regression.cpp)
    target_link_libraries(test_kinematics_regression PRIVATE rk4)
    target_link_libraries(test_detector_regression PRIVATE rk4)
endif()
if(BUILD_TESTING)
    add_executable(test_stopping_cli test/stopping_cli.cpp)
    target_link_libraries(test_stopping_cli PRIVATE rk4)
endif()
if(BUILD_TESTING)
    add_executable(test_h5_reader_cli test/h5_reader_cli.cpp)
    target_link_libraries(test_h5_reader_cli PRIVATE rk4)
endif()
install(FILES notices/numpy-random.txt DESTINATION ${CMAKE_INSTALL_DATADIR}/spyralpp/licenses)
install(FILES notices/spyral-reference.txt DESTINATION ${CMAKE_INSTALL_DATADIR}/spyralpp/licenses)
if(BUILD_TESTING)
    add_executable(test_pcg64_cli test/pcg64_cli.cpp)
    target_include_directories(test_pcg64_cli PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
endif()
if(BUILD_TESTING)
    add_executable(test_waveform_cli test/waveform_cli.cpp)
    target_link_libraries(test_waveform_cli PRIVATE rk4)
endif()
if(BUILD_TESTING)
    add_executable(test_hdbscan_cli test/hdbscan_cli.cpp)
    target_link_libraries(test_hdbscan_cli PRIVATE rk4)
endif()
if(BUILD_TESTING)
    add_executable(test_clusterer_snapshot_regression test/clusterer_snapshot_regression.cpp)
    target_link_libraries(test_clusterer_snapshot_regression PRIVATE rk4)
    add_test(NAME clusterer_snapshot_regression COMMAND test_clusterer_snapshot_regression)
endif()
install(FILES notices/scikit-learn-hdbscan.txt DESTINATION ${CMAKE_INSTALL_DATADIR}/spyralpp/licenses)

# Pin NumPy's SIMD argsort kernel: equal-distance MST ties affect HDBSCAN labels.
target_sources(rk4 PRIVATE src/reco/NumpyArgsort.C)
if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64" AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_sources(rk4 PRIVATE src/reco/NumpyArgsortAvx2.C src/reco/NumpyArgsortAvx512.C)
    set_source_files_properties(src/reco/NumpyArgsortAvx2.C PROPERTIES COMPILE_OPTIONS "-mavx2")
    set_source_files_properties(src/reco/NumpyArgsortAvx512.C PROPERTIES COMPILE_OPTIONS "-mavx2;-mavx512f;-mavx512dq;-mavx512vl")
    target_compile_definitions(rk4 PRIVATE RK4_NUMPY_X86_ARGSORT)
    target_include_directories(rk4 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/extern/x86-simd-sort/src)
endif()
install(FILES extern/x86-simd-sort/LICENSE.md DESTINATION ${CMAKE_INSTALL_DATADIR}/spyralpp/licenses RENAME x86-simd-sort-LICENSE.md)
if(BUILD_TESTING)
    add_executable(test_smoothing_spline_cli test/smoothing_spline_cli.cpp)
    target_include_directories(test_smoothing_spline_cli PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
endif()

if(BUILD_TESTING)
    add_executable(test_reco_config test/testRecoConfig.cpp)
    target_link_libraries(test_reco_config PRIVATE rk4)
    add_test(NAME reco_config COMMAND test_reco_config)
endif()

if(BUILD_TESTING)
    add_executable(test_pid_state test/testPidState.cpp)
    target_include_directories(test_pid_state PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
    add_test(NAME pid_state COMMAND test_pid_state)
endif()

if(BUILD_TESTING)
    add_executable(test_result_io_state_metadata test/testResultIOStateMetadata.cpp)
    target_link_libraries(test_result_io_state_metadata PRIVATE rk4)
    add_test(NAME result_io_state_metadata COMMAND test_result_io_state_metadata)
endif()

if(BUILD_TESTING)
    add_executable(test_workspace_io test/testWorkspaceIO.cpp)
    target_link_libraries(test_workspace_io PRIVATE rk4)
    add_test(NAME workspace_io COMMAND test_workspace_io)
    set_tests_properties(workspace_io PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif()
if(BUILD_TESTING AND SPYRALPP_BUILD_APPS)
    add_executable(test_workspace_solver test/testWorkspaceSolver.cpp)
    target_link_libraries(test_workspace_solver PRIVATE rk4)
    add_test(NAME workspace_solver COMMAND test_workspace_solver $<TARGET_FILE:rk4_workspace>)
    set_tests_properties(workspace_solver PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif()
