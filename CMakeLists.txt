cmake_minimum_required(VERSION 3.30)

# --- vcpkg tool pin (ARCHITECTURE.md, Q3) ---------------------------------------
# The vcpkg tool checkout named by VCPKG_ROOT must be at the same commit as the
# manifest's builtin-baseline, so the tool and the port versions cannot drift apart.
if(NOT CMAKE_HOST_WIN32)
  message(FATAL_ERROR "NiTCAD is Windows-only.")
endif()
if(NOT DEFINED ENV{VCPKG_ROOT})
  message(FATAL_ERROR "VCPKG_ROOT is not set. Point it at the pinned standalone vcpkg checkout.")
endif()
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/vcpkg.json" NITCAD_MANIFEST)
string(JSON NITCAD_VCPKG_COMMIT GET "${NITCAD_MANIFEST}" "builtin-baseline")
execute_process(
  COMMAND git -C "$ENV{VCPKG_ROOT}" rev-parse HEAD
  OUTPUT_VARIABLE NITCAD_VCPKG_HEAD
  RESULT_VARIABLE NITCAD_VCPKG_RC
  OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
if(NOT NITCAD_VCPKG_RC EQUAL 0 OR NOT NITCAD_VCPKG_HEAD STREQUAL NITCAD_VCPKG_COMMIT)
  message(FATAL_ERROR
    "VCPKG_ROOT ($ENV{VCPKG_ROOT}) is not a git checkout at the pinned commit "
    "${NITCAD_VCPKG_COMMIT} (found: '${NITCAD_VCPKG_HEAD}').")
endif()

# Dynamic CRT, matching the x64-windows-static-md triplet (N1).
set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")

project(NiTCAD VERSION 0.0.0 LANGUAGES CXX)

# --- Toolchain gate (ARCHITECTURE.md, D1 and Q6) ----------------------------------
# MSVC only. 14.51 (compiler 19.51) is the verified reference toolset.
if(NOT MSVC)
  message(FATAL_ERROR "NiTCAD requires MSVC (found ${CMAKE_CXX_COMPILER_ID}).")
endif()
if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 19.51)
  message(FATAL_ERROR
    "MSVC ${CMAKE_CXX_COMPILER_VERSION} is older than the verified reference toolset 19.51 (14.51).")
endif()

# The language mode is passed explicitly. CMAKE_CXX_STANDARD is deliberately not used:
# CMake maps 23 to /std:c++latest, a rolling mode that is not C++23.
# /options:strict turns an unknown option (for example an ignored /std:c++23, which
# otherwise falls back to C++14 with only warning D9002) into an error.
set(NITCAD_LANG_FLAGS /std:c++23preview /options:strict /EHsc /utf-8 /permissive- /fp:precise)

# Configure-time probe: the compiler must really be in C++23 mode and provide the
# library features the code uses.
unset(NITCAD_CXX23_GATE CACHE)
file(WRITE "${CMAKE_BINARY_DIR}/nitcad_cxx23_gate.cpp" [=[
#include <version>
#include <expected>
#include <stop_token>
#if _MSVC_LANG != 202302L
#error "NiTCAD requires /std:c++23preview (_MSVC_LANG == 202302L)"
#endif
#ifndef __cpp_lib_expected
#error "__cpp_lib_expected is not defined"
#endif
#ifndef __cpp_lib_jthread
#error "__cpp_lib_jthread is not defined"
#endif
int main() { return 0; }
]=])
try_compile(NITCAD_CXX23_GATE
  SOURCES "${CMAKE_BINARY_DIR}/nitcad_cxx23_gate.cpp"
  COMPILE_DEFINITIONS ${NITCAD_LANG_FLAGS}
  OUTPUT_VARIABLE NITCAD_CXX23_GATE_OUTPUT)
if(NOT NITCAD_CXX23_GATE)
  message(FATAL_ERROR "C++23 gate failed:\n${NITCAD_CXX23_GATE_OUTPUT}")
endif()
message(STATUS "C++23 gate passed: /std:c++23preview, MSVC ${CMAKE_CXX_COMPILER_VERSION}")

# Warning policy: first-party code at /W4 with warnings as errors; imported headers
# are system headers and their warnings are silenced. CMake's default /W3 is removed so
# it does not conflict with /W4.
string(REGEX REPLACE "/W[0-4]" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
add_compile_options(${NITCAD_LANG_FLAGS} /W4 /WX /external:W0)
add_compile_definitions(UNICODE _UNICODE NOMINMAX WIN32_LEAN_AND_MEAN)

# --- Dependencies (from vcpkg.json) ------------------------------------------------
find_package(Catch2 3 CONFIG REQUIRED)
find_package(Eigen3 CONFIG REQUIRED)
# Optional: the build stamps its git revision into a run's provenance (Unit 42); "unknown" without.
find_package(Git QUIET)

# --- Scaffold test -----------------------------------------------------------------
enable_testing()
add_executable(nitcad_scaffold_test tests/scaffold_test.cpp)
target_link_libraries(nitcad_scaffold_test PRIVATE Catch2::Catch2WithMain)
add_test(NAME scaffold COMMAND nitcad_scaffold_test)

# --- base layer (constants, Error, NITCAD_EXPECTS) ---------------------------------
add_library(nitcad_base STATIC src/base/contract.cpp)
add_library(NiTCAD::base ALIAS nitcad_base)
target_include_directories(nitcad_base PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")

add_executable(nitcad_base_test
  tests/base/constants_test.cpp
  tests/base/error_test.cpp
  tests/base/contract_test.cpp)
target_link_libraries(nitcad_base_test PRIVATE NiTCAD::base Catch2::Catch2WithMain)
add_test(NAME base COMMAND nitcad_base_test)

# A violated NITCAD_EXPECTS terminates the process, so contract_test.cpp runs this probe
# executable as a child process and inspects its exit code and stderr.
add_executable(nitcad_contract_probe tests/base/contract_probe.cpp)
target_link_libraries(nitcad_contract_probe PRIVATE NiTCAD::base)
add_dependencies(nitcad_base_test nitcad_contract_probe)
target_compile_definitions(nitcad_base_test PRIVATE
  NITCAD_CONTRACT_PROBE_PATH="$<TARGET_FILE:nitcad_contract_probe>")

# --- MKL runtime for the PARDISO tests (Unit 18) ------------------------------------
# MKL is optional and never linked: the mkl_pardiso backend loads mkl_rt at run time from a path
# the program gives (linalg/mkl_runtime.hpp). The tests take that path from this cache variable;
# when it is empty they skip the PARDISO cases, unless NITCAD_REQUIRE_MKL is on (CI), which makes
# a missing runtime a configuration error.
set(NITCAD_MKL_RUNTIME "" CACHE FILEPATH "mkl_rt DLL for the PARDISO tests (optional)")
option(NITCAD_REQUIRE_MKL "Fail the configuration when NITCAD_MKL_RUNTIME is not a file" OFF)
if(NITCAD_REQUIRE_MKL AND NOT EXISTS "${NITCAD_MKL_RUNTIME}")
  message(FATAL_ERROR "NITCAD_REQUIRE_MKL is on but NITCAD_MKL_RUNTIME ('${NITCAD_MKL_RUNTIME}') is not a file.")
endif()
file(TO_CMAKE_PATH "${NITCAD_MKL_RUNTIME}" NITCAD_MKL_RUNTIME_PATH)

# --- linalg layer (sparse matrix, backend-neutral solver, Eigen SparseLU, PARDISO and GMRES-ILUT backends) ---
add_library(nitcad_linalg STATIC
  src/linalg/sparse_matrix.cpp
  src/linalg/linear_solver.cpp
  src/linalg/eigen_sparse_lu.cpp
  src/linalg/eigen_ilut.cpp
  src/linalg/mkl_runtime.cpp
  src/linalg/mkl_pardiso.cpp)
add_library(NiTCAD::linalg ALIAS nitcad_linalg)
target_include_directories(nitcad_linalg PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
# Eigen is PRIVATE so its include path never reaches linalg consumers (6.10, Q5);
# tests/linalg/header_boundary_test.cpp fails to compile if it does.
target_link_libraries(nitcad_linalg PUBLIC NiTCAD::base PRIVATE Eigen3::Eigen)

add_executable(nitcad_linalg_test
  tests/linalg/sparse_matrix_test.cpp
  tests/linalg/linear_solver_test.cpp
  tests/linalg/pardiso_test.cpp
  tests/linalg/gmres_kernel_test.cpp
  tests/linalg/gmres_ilut_test.cpp
  tests/linalg/header_boundary_test.cpp
  tests/linalg/compensated_dot_test.cpp)
target_link_libraries(nitcad_linalg_test PRIVATE NiTCAD::linalg Catch2::Catch2WithMain)
target_compile_definitions(nitcad_linalg_test PRIVATE
  NITCAD_MKL_RUNTIME_PATH="${NITCAD_MKL_RUNTIME_PATH}")
# A process that never loads MKL: mkl_pardiso must refuse to start (its own process, because the
# PARDISO tests load the runtime and Catch2 runs tests in random order).
add_executable(nitcad_mkl_absent_probe tests/linalg/mkl_absent_probe.cpp)
target_link_libraries(nitcad_mkl_absent_probe PRIVATE NiTCAD::linalg)
add_test(NAME linalg_mkl_absent COMMAND nitcad_mkl_absent_probe)
# One process per logical CPU, each pinned: the Eigen backend gives the same bits on every core
# (Unit 41; its own executable, because Eigen detects the cache sizes once per process).
add_executable(nitcad_eigen_blocking_probe tests/linalg/eigen_blocking_probe.cpp)
target_link_libraries(nitcad_eigen_blocking_probe PRIVATE NiTCAD::linalg)
add_test(NAME linalg_eigen_blocking COMMAND nitcad_eigen_blocking_probe)
add_test(NAME linalg COMMAND nitcad_linalg_test)

# --- mesh layer (node/edge/control-volume graph, tensor-grid producer) -------------
add_library(nitcad_mesh STATIC
  src/mesh/mesh.cpp
  src/mesh/tensor_grid.cpp
  src/mesh/tensor_cells.cpp
  src/mesh/exact_integer.cpp
  src/mesh/predicates.cpp
  src/mesh/dual_geometry.cpp
  src/mesh/simplex_mesh.cpp
  src/mesh/interfaces.cpp
  src/mesh/triangulation2d.cpp
  src/mesh/triangulation3d.cpp
  src/mesh/generate2d.cpp
  src/mesh/generate3d.cpp)
add_library(NiTCAD::mesh ALIAS nitcad_mesh)
target_include_directories(nitcad_mesh PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(nitcad_mesh PUBLIC NiTCAD::base)

add_executable(nitcad_mesh_test
  tests/mesh/mesh_test.cpp
  tests/mesh/tensor_grid_test.cpp
  tests/mesh/tensor_cells_test.cpp
  tests/mesh/predicates_test.cpp
  tests/mesh/filter_test.cpp
  tests/mesh/simplex_mesh_test.cpp
  tests/mesh/interfaces_test.cpp
  tests/mesh/origin_test.cpp
  tests/mesh/generate2d_test.cpp
  tests/mesh/triangulation3d_test.cpp
  tests/mesh/generate3d_test.cpp)
# The predicate and dual-geometry tests read the internal headers directly (Unit 16).
target_include_directories(nitcad_mesh_test PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src/mesh")
target_link_libraries(nitcad_mesh_test PRIVATE NiTCAD::mesh Catch2::Catch2WithMain)
add_test(NAME mesh COMMAND nitcad_mesh_test)

# --- physics layer (silicon parameters, Caughey-Thomas, SRH, Boltzmann statistics) --
add_library(nitcad_physics STATIC
  src/physics/semiconductor.cpp
  src/physics/mobility.cpp
  src/physics/recombination.cpp
  src/physics/statistics.cpp
  src/physics/bandgap_narrowing.cpp
  src/physics/field_mobility.cpp
  src/physics/fermi_dirac.cpp
  src/physics/thermionic_emission.cpp
  src/physics/ionization.cpp
  src/physics/insulator.cpp
  src/physics/interface_traps.cpp
  src/physics/impact_ionization.cpp
  src/physics/band_to_band.cpp
  src/physics/thermal.cpp)
add_library(NiTCAD::physics ALIAS nitcad_physics)
target_include_directories(nitcad_physics PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(nitcad_physics PUBLIC NiTCAD::base)

add_executable(nitcad_physics_test
  tests/physics/semiconductor_test.cpp
  tests/physics/mobility_test.cpp
  tests/physics/recombination_test.cpp
  tests/physics/statistics_test.cpp
  tests/physics/heavy_doping_test.cpp
  tests/physics/field_mobility_test.cpp
  tests/physics/fermi_dirac_test.cpp
  tests/physics/materials_test.cpp
  tests/physics/interface_traps_test.cpp
  tests/physics/impact_ionization_test.cpp
  tests/physics/band_to_band_test.cpp
  tests/physics/thermal_test.cpp
  tests/physics/masetti_test.cpp
  tests/physics/band_set_test.cpp
  tests/physics/ionization_density_test.cpp)
target_link_libraries(nitcad_physics_test PRIVATE NiTCAD::physics Catch2::Catch2WithMain)
add_test(NAME physics COMMAND nitcad_physics_test)

# --- device layer (regions, doping, contacts; validated problem description) -------
add_library(nitcad_device STATIC
  src/device/device.cpp
  src/device/doping_profiles.cpp)
add_library(NiTCAD::device ALIAS nitcad_device)
target_include_directories(nitcad_device PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(nitcad_device PUBLIC NiTCAD::base NiTCAD::mesh NiTCAD::physics)

add_executable(nitcad_device_test
  tests/device/device_test.cpp
  tests/device/topology_test.cpp
  tests/device/insulator_test.cpp
  tests/device/certificate_test.cpp
  tests/device/doping_profiles_test.cpp)
target_link_libraries(nitcad_device_test PRIVATE NiTCAD::device Catch2::Catch2WithMain)
add_test(NAME device COMMAND nitcad_device_test)

# --- assemble layer (scaling, Bernoulli/SG flux, ohmic and gate terms, Poisson, DD) ---
add_library(nitcad_assemble STATIC
  src/assemble/scaling.cpp
  src/assemble/scaled_device.cpp
  src/assemble/gate.cpp
  src/assemble/interface_edges.cpp
  src/assemble/contact_bias.cpp
  src/assemble/equilibrium_poisson.cpp
  src/assemble/drift_diffusion.cpp
  src/assemble/electrothermal.cpp
  src/assemble/tunnel_paths.cpp
  src/assemble/observer.cpp
  src/assemble/fields.cpp
  src/assemble/terms.cpp
  src/assemble/thermal_terms.cpp
  src/assemble/composition.cpp)
add_library(NiTCAD::assemble ALIAS nitcad_assemble)
target_include_directories(nitcad_assemble PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(nitcad_assemble PUBLIC
  NiTCAD::base NiTCAD::linalg NiTCAD::mesh NiTCAD::physics NiTCAD::device)

add_executable(nitcad_assemble_test
  tests/assemble/flux_test.cpp
  tests/assemble/poisson_test.cpp
  tests/assemble/drift_diffusion_test.cpp
  tests/assemble/gate_test.cpp
  tests/assemble/heterojunction_test.cpp
  tests/assemble/insulator_test.cpp
  tests/assemble/transient_test.cpp
  tests/assemble/small_signal_test.cpp
  tests/assemble/impact_ionization_test.cpp
  tests/assemble/band_to_band_test.cpp
  tests/assemble/electrothermal_test.cpp
  tests/assemble/fields_test.cpp
  tests/assemble/terms_test.cpp
  tests/baseline/trace_listener.cpp)
target_link_libraries(nitcad_assemble_test PRIVATE NiTCAD::assemble Catch2::Catch2WithMain)
# The terms (Unit 40) are internal to the assemble layer; their own test reads them.
target_include_directories(nitcad_assemble_test PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src/assemble")
add_test(NAME assemble COMMAND nitcad_assemble_test)

# --- results layer (plain result data; header-only) --------------------------------
add_library(nitcad_results INTERFACE)
add_library(NiTCAD::results ALIAS nitcad_results)
target_include_directories(nitcad_results INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(nitcad_results INTERFACE NiTCAD::base)

# --- analysis layer (extraction from results; never solves) ------------------------
add_library(nitcad_analysis STATIC
  src/analysis/curve.cpp
  src/analysis/dc.cpp
  src/analysis/cv.cpp
  src/analysis/two_port.cpp)
add_library(NiTCAD::analysis ALIAS nitcad_analysis)
target_include_directories(nitcad_analysis PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(nitcad_analysis PUBLIC NiTCAD::base NiTCAD::results)

add_executable(nitcad_analysis_test
  tests/analysis/curve_test.cpp
  tests/analysis/dc_test.cpp
  tests/analysis/cv_test.cpp
  tests/analysis/two_port_test.cpp)
target_link_libraries(nitcad_analysis_test PRIVATE NiTCAD::analysis Catch2::Catch2WithMain)
add_test(NAME analysis COMMAND nitcad_analysis_test)

# --- io layer (files: meshes, sweeps, checkpoints; Unit 29, approved 2026-10-09) ----------
add_library(nitcad_io STATIC
  src/io/record.cpp
  src/io/files.cpp
  src/io/device_file.cpp)
add_library(NiTCAD::io ALIAS nitcad_io)
target_include_directories(nitcad_io PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(nitcad_io PUBLIC NiTCAD::base NiTCAD::mesh NiTCAD::device NiTCAD::results)

add_executable(nitcad_io_test
  tests/io/record_test.cpp
  tests/io/files_test.cpp
  tests/io/origin_test.cpp
  tests/io/generated_test.cpp
  tests/io/device_file_test.cpp
  tests/io/provenance_test.cpp
  tests/io/run_files_test.cpp
  tests/io/rates_files_test.cpp)
# The record test reads the private encoding directly.
target_include_directories(nitcad_io_test PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src/io")
target_link_libraries(nitcad_io_test PRIVATE NiTCAD::io Catch2::Catch2WithMain)
# Files written by earlier builds, which must keep loading (Unit 42, PG3).
target_compile_definitions(nitcad_io_test PRIVATE
  NITCAD_IO_FIXTURES="${CMAKE_CURRENT_SOURCE_DIR}/tests/io/fixtures")
add_test(NAME io COMMAND nitcad_io_test)

# --- render layer (the views: CPU scene, then drawing; Unit 25a, DECISIONS.md VR1-VR9, VQ8) ------
# HLSL compiled at build time with the SDK's dxc into headers of bytecode (VQ2).
find_program(NITCAD_DXC dxc
  HINTS "$ENV{WindowsSdkVerBinPath}/x64" "C:/Program Files (x86)/Windows Kits/10/bin/10.0.26100.0/x64"
  REQUIRED)
set(NITCAD_SHADER_DIR "${CMAKE_CURRENT_BINARY_DIR}/shaders")
set(NITCAD_SHADER_HEADERS "")
foreach(entry
    "scene mesh_vs vs_6_0" "scene opaque_ps ps_6_0" "scene transparent_ps ps_6_0"
    "scene line_vs vs_6_0" "scene line_ps ps_6_0"
    "post fullscreen_vs vs_6_0" "post composite_ps ps_6_0" "post resolve_ps ps_6_0")
  separate_arguments(parts NATIVE_COMMAND "${entry}")
  list(GET parts 0 file)
  list(GET parts 1 name)
  list(GET parts 2 profile)
  set(header "${NITCAD_SHADER_DIR}/${name}.h")
  add_custom_command(OUTPUT "${header}"
    COMMAND ${CMAKE_COMMAND} -E make_directory "${NITCAD_SHADER_DIR}"
    COMMAND "${NITCAD_DXC}" -nologo -WX -O3 -T ${profile} -E ${name} -Vn ${name}_bytecode
            -Fh "${header}" "${CMAKE_CURRENT_SOURCE_DIR}/src/render/shaders/${file}.hlsl"
    DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/render/shaders/${file}.hlsl"
    COMMENT "dxc ${file}.hlsl ${name}")
  list(APPEND NITCAD_SHADER_HEADERS "${header}")
endforeach()
add_library(nitcad_render STATIC
  src/render/color.cpp
  src/render/geometry.cpp
  src/render/scene.cpp
  src/render/vectors.cpp
  src/render/streamlines.cpp
  src/render/cut_line.cpp
  src/render/quality.cpp
  src/render/compare.cpp
  src/render/renderer.cpp
  src/render/png.cpp
  ${NITCAD_SHADER_HEADERS})
add_library(NiTCAD::render ALIAS nitcad_render)
target_include_directories(nitcad_render PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(nitcad_render PUBLIC NiTCAD::base NiTCAD::mesh NiTCAD::results
  NiTCAD::analysis NiTCAD::device PRIVATE d3d12 dxgi d3d11 d2d1 dwrite windowscodecs ole32)
target_include_directories(nitcad_render PRIVATE "${NITCAD_SHADER_DIR}")

add_executable(nitcad_render_test
  tests/render/color_test.cpp
  tests/render/geometry_test.cpp
  tests/render/scene_test.cpp
  tests/render/vectors_test.cpp
  tests/render/iso_test.cpp
  tests/render/streamlines_test.cpp
  tests/render/arrows_test.cpp
  tests/render/clip_test.cpp
  tests/render/exploded_test.cpp
  tests/render/film_test.cpp
  tests/render/cut_line_test.cpp
  tests/render/quality_test.cpp
  tests/render/compare_test.cpp
  tests/render/image_test.cpp
  tests/render/window_test.cpp)
# The image test solves its MOSFETs (tests may use any layer; render itself may not).
target_link_libraries(nitcad_render_test PRIVATE NiTCAD::render NiTCAD::solve Catch2::Catch2WithMain)
target_compile_definitions(nitcad_render_test PRIVATE
  NITCAD_RENDER_REFERENCES="${CMAKE_CURRENT_SOURCE_DIR}/tests/render/references")
add_test(NAME render COMMAND nitcad_render_test)
# VG1: render's includes stay within base, mesh, results, analysis and device.
add_test(NAME render_layering
  COMMAND ${CMAKE_COMMAND} -DROOT=${CMAKE_CURRENT_SOURCE_DIR} -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/render/layering.cmake)

# --- solve layer (Newton, equilibrium and bias solves, sweeps, cancellation, progress) ---
add_library(nitcad_solve STATIC
  src/solve/equilibrium.cpp
  src/solve/bias.cpp
  src/solve/run_record.cpp
  src/solve/weighted_current.cpp
  src/solve/waveform.cpp
  src/solve/transient.cpp
  src/solve/small_signal.cpp
  src/solve/trace.cpp
  src/solve/provenance.cpp)
add_library(NiTCAD::solve ALIAS nitcad_solve)
target_include_directories(nitcad_solve PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(nitcad_solve PUBLIC
  NiTCAD::base NiTCAD::linalg NiTCAD::assemble NiTCAD::device NiTCAD::results)
# The build facts of a run's provenance (Unit 42; DECISIONS.md PV3), rewritten at every build by
# cmake/build_info.cmake but only when they change, so only provenance.cpp ever recompiles.
set(NITCAD_BUILD_INFO "${CMAKE_CURRENT_BINARY_DIR}/generated/NiTCAD/build_info.hpp")
add_custom_target(nitcad_build_info
  COMMAND "${CMAKE_COMMAND}" "-DSOURCE=${CMAKE_CURRENT_SOURCE_DIR}" "-DOUTPUT=${NITCAD_BUILD_INFO}"
    "-DGIT=${GIT_EXECUTABLE}" "-DVERSION=${PROJECT_VERSION}" "-DBUILD_TYPE=${CMAKE_BUILD_TYPE}"
    "-DCOMPILER=MSVC ${CMAKE_CXX_COMPILER_VERSION}"
    -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/build_info.cmake"
  BYPRODUCTS "${NITCAD_BUILD_INFO}"
  VERBATIM)
add_dependencies(nitcad_solve nitcad_build_info)
target_include_directories(nitcad_solve PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/generated")

add_executable(nitcad_solve_test
  tests/solve/newton_test.cpp
  tests/solve/equilibrium_test.cpp
  tests/solve/bias_test.cpp
  tests/solve/step_control_test.cpp
  tests/solve/sweep_test.cpp
  tests/solve/heavy_doping_test.cpp
  tests/solve/mos_test.cpp
  tests/solve/mosfet_test.cpp
  tests/solve/unstructured_test.cpp
  tests/solve/unstructured3d_test.cpp
  tests/solve/field_mobility_test.cpp
  tests/solve/fermi_dirac_test.cpp
  tests/solve/heterojunction_test.cpp
  tests/solve/insulator_test.cpp
  tests/solve/transient_test.cpp
  tests/solve/small_signal_test.cpp
  tests/solve/impact_ionization_test.cpp
  tests/solve/btbt_test.cpp
  tests/solve/electrothermal_test.cpp
  tests/solve/restart_test.cpp
  tests/solve/weighted_current_test.cpp
  tests/solve/hole_mobility_test.cpp
  tests/solve/default_physics_test.cpp
  tests/solve/provenance_test.cpp
  tests/solve/rates_test.cpp
  tests/solve/process_memory.cpp
  tests/baseline/trace_listener.cpp)
# The solve tests also read their results through the analysis extractors (Unit 24 gates).
target_link_libraries(nitcad_solve_test PRIVATE NiTCAD::solve NiTCAD::analysis NiTCAD::io
  Catch2::Catch2WithMain)
# The 3D sliver stress test meshes perturbed lattices with the generator's Delaunay (Unit 16).
target_include_directories(nitcad_solve_test PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src/mesh")
target_compile_definitions(nitcad_solve_test PRIVATE
  NITCAD_MKL_RUNTIME_PATH="${NITCAD_MKL_RUNTIME_PATH}"
  NITCAD_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}" NITCAD_GIT="${GIT_EXECUTABLE}")
add_test(NAME solve COMMAND nitcad_solve_test)
# The 2D MOSFET gates (tests/solve/mosfet_test.cpp, hidden tag [.mosfet]) take about a minute in
# Release and far longer in Debug, so they run in the Release configuration only.
if(CMAKE_BUILD_TYPE STREQUAL "Release")
  add_test(NAME solve_mosfet COMMAND nitcad_solve_test "[mosfet]")
  # Likewise the slowest transient gates (tests/solve/transient_test.cpp, hidden tag [.transient]).
  add_test(NAME solve_transient COMMAND nitcad_solve_test "[transient]")
  # And the breakdown gates (tests/solve/impact_ionization_test.cpp, hidden tag [.breakdown]).
  add_test(NAME solve_breakdown COMMAND nitcad_solve_test "[breakdown]")
  # And the slow band-to-band gates (tests/solve/btbt_test.cpp, hidden tag [.btbt]).
  add_test(NAME solve_btbt COMMAND nitcad_solve_test "[btbt]")
  # And the device gates on generated 2D meshes (tests/solve/unstructured_test.cpp, [.unstructured]).
  add_test(NAME solve_unstructured COMMAND nitcad_solve_test "[unstructured]")
endif()

# --- deck layer (the text input deck, its meaning and the batch runner; Unit 31) ----------
add_library(nitcad_deck STATIC
  src/deck/syntax.cpp
  src/deck/build.cpp
  src/deck/run.cpp)
add_library(NiTCAD::deck ALIAS nitcad_deck)
target_include_directories(nitcad_deck PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(nitcad_deck PUBLIC NiTCAD::base NiTCAD::io NiTCAD::solve)

add_executable(nitcad_deck_test
  tests/deck/syntax_test.cpp
  tests/deck/build_test.cpp
  tests/deck/run_test.cpp
  tests/deck/analyses_test.cpp
  tests/deck/analyses_run_test.cpp
  tests/baseline/trace_listener.cpp)
target_include_directories(nitcad_deck_test PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src/deck")
# The plot tests read the images back (Unit 25a).
target_link_libraries(nitcad_deck_test PRIVATE NiTCAD::deck NiTCAD::render Catch2::Catch2WithMain)
# The run tests start the nitcad executable for the exit status and Ctrl+C-free paths.
add_dependencies(nitcad_deck_test nitcad)
target_compile_definitions(nitcad_deck_test PRIVATE
  NITCAD_CLI_PATH="$<TARGET_FILE:nitcad>"
  NITCAD_IO_FIXTURES="${CMAKE_CURRENT_SOURCE_DIR}/tests/io/fixtures"
  NITCAD_MKL_RUNTIME_PATH="${NITCAD_MKL_RUNTIME_PATH}")
add_test(NAME deck COMMAND nitcad_deck_test)

# The batch runner (`nitcad check|run <deck>`).
add_executable(nitcad src/cli/main.cpp src/cli/plot.cpp)
target_link_libraries(nitcad PRIVATE NiTCAD::deck NiTCAD::render)
# `nitcad --version` names the build (Unit 42, PV6): a 40-digit revision, or unknown without git.
add_test(NAME cli_version COMMAND nitcad --version)
set_tests_properties(cli_version PROPERTIES PASS_REGULAR_EXPRESSION
  "^nitcad ${PROJECT_VERSION} \\(revision ([0-9a-f]+|unknown)(, modified)?, ${CMAKE_BUILD_TYPE}, MSVC [0-9.]+, Eigen [0-9]")

# --- app layer (the Windows application; Unit 26, DECISIONS.md VA1-VA9, AQ1-AQ4, WG1-WG8) ---------
# The logic, plain C++ under test (26-1): documents, the run controller, the workspace, view settings,
# the menus' deck sections, the per-user state.
add_library(nitcad_app STATIC
  src/app/document.cpp
  src/app/run_controller.cpp
  src/app/workspace.cpp
  src/app/view_settings.cpp
  src/app/deck_sections.cpp
  src/app/app_state.cpp
  src/app/commands.cpp
  src/app/layout.cpp
  src/app/probe.cpp
  src/app/examples.cpp
  src/app/run_monitor.cpp
  src/app/deck_vocabulary.cpp
  src/app/live_check.cpp
  src/app/compare.cpp)
add_library(NiTCAD::app ALIAS nitcad_app)
target_include_directories(nitcad_app PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(nitcad_app PUBLIC NiTCAD::deck NiTCAD::render NiTCAD::analysis NiTCAD::io)

add_executable(nitcad_app_test
  tests/app/document_test.cpp
  tests/app/run_controller_test.cpp
  tests/app/workspace_test.cpp
  tests/app/view_settings_test.cpp
  tests/app/deck_sections_test.cpp
  tests/app/app_state_test.cpp
  tests/app/commands_test.cpp
  tests/app/layout_test.cpp
  tests/app/probe_test.cpp
  tests/app/examples_test.cpp
  tests/app/run_monitor_test.cpp
  tests/app/deck_vocabulary_test.cpp
  tests/app/live_check_test.cpp
  tests/app/compare_test.cpp)
target_link_libraries(nitcad_app_test PRIVATE NiTCAD::app Catch2::Catch2WithMain user32)
# The view tests compare the app's images with nitcad plot's (WG4's choices).
add_dependencies(nitcad_app_test nitcad)
target_compile_definitions(nitcad_app_test PRIVATE NITCAD_CLI_PATH="$<TARGET_FILE:nitcad>"
  NITCAD_EXAMPLES_DIR="${CMAKE_CURRENT_SOURCE_DIR}/examples"
  NITCAD_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
add_test(NAME app COMMAND nitcad_app_test)
# The window (26-2): the main window, the view and the deck pane on Win32, in a library the window
# tests share with the executable.
add_library(nitcad_app_window STATIC
  src/app/win32.cpp
  src/app/viewport.cpp
  src/app/deck_pane.cpp
  src/app/main_window.cpp
  src/app/pages.cpp
  src/app/self_test.cpp)
target_link_libraries(nitcad_app_window PUBLIC NiTCAD::app comctl32 comdlg32 shell32 ole32 user32 gdi32 dxgi dxguid)
# The start page lists the examples beside the executable, else the source tree's (26-3, VQ11).
target_compile_definitions(nitcad_app_window PRIVATE NITCAD_EXAMPLES_DIR="${CMAKE_CURRENT_SOURCE_DIR}/examples")
target_include_directories(nitcad_app_window PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src/app")
# nitcad-app (AQ4): a separate GUI executable; the nitcad batch runner keeps no window code.
add_executable(nitcad-app WIN32 app/main.cpp)
target_link_libraries(nitcad-app PRIVATE nitcad_app_window)
add_executable(nitcad_app_window_test tests/app/window_test.cpp tests/app/self_test_test.cpp
  tests/app/analysis_window_test.cpp)
target_link_libraries(nitcad_app_window_test PRIVATE nitcad_app_window Catch2::Catch2WithMain dxgi dxguid)
target_compile_definitions(nitcad_app_window_test PRIVATE NITCAD_CLI_PATH="$<TARGET_FILE:nitcad>"
  NITCAD_APP_PATH="$<TARGET_FILE:nitcad-app>" NITCAD_EXAMPLES_DIR="${CMAKE_CURRENT_SOURCE_DIR}/examples")
# The self-test (WG4) runs nitcad-app and nitcad plot.
add_dependencies(nitcad_app_window_test nitcad-app nitcad)
add_test(NAME app_window COMMAND nitcad_app_window_test)
# The 3D MOSFET example's run and self-tests (WG8, WG4) take minutes in Debug: Release only, as the other
# slow gates (hidden tags [.examples3d] and [.selftest3d]).
if(CMAKE_BUILD_TYPE STREQUAL "Release")
  add_test(NAME app_examples_3d COMMAND nitcad_app_test "[examples3d]")
  add_test(NAME app_window_3d COMMAND nitcad_app_window_test "[selftest3d]")
endif()
# WG6: app includes every layer but physics and assemble (ARCHITECTURE.md 3).
add_test(NAME app_layering
  COMMAND ${CMAKE_COMMAND} -DROOT=${CMAKE_CURRENT_SOURCE_DIR} -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/app/layering.cmake)

# --- validation against published data (Unit 33; benchmarks/ decks) ----------------------
add_executable(nitcad_validation_test tests/validation/published_test.cpp
  tests/baseline/trace_listener.cpp)
target_link_libraries(nitcad_validation_test PRIVATE NiTCAD::deck Catch2::Catch2WithMain)
target_compile_definitions(nitcad_validation_test PRIVATE
  NITCAD_BENCHMARKS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/benchmarks")
# The benchmarks solve many devices and pin measured gaps; they run in the Release configuration.
if(CMAKE_BUILD_TYPE STREQUAL "Release")
  add_test(NAME validation COMMAND nitcad_validation_test)
endif()
