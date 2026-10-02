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

# --- linalg layer (sparse matrix, backend-neutral solver, Eigen SparseLU backend) ---
add_library(nitcad_linalg STATIC
  src/linalg/sparse_matrix.cpp
  src/linalg/linear_solver.cpp
  src/linalg/eigen_sparse_lu.cpp)
add_library(NiTCAD::linalg ALIAS nitcad_linalg)
target_include_directories(nitcad_linalg PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
# Eigen is PRIVATE so its include path never reaches linalg consumers (6.10, Q5);
# tests/linalg/header_boundary_test.cpp fails to compile if it does.
target_link_libraries(nitcad_linalg PUBLIC NiTCAD::base PRIVATE Eigen3::Eigen)

add_executable(nitcad_linalg_test
  tests/linalg/sparse_matrix_test.cpp
  tests/linalg/linear_solver_test.cpp
  tests/linalg/header_boundary_test.cpp)
target_link_libraries(nitcad_linalg_test PRIVATE NiTCAD::linalg Catch2::Catch2WithMain)
add_test(NAME linalg COMMAND nitcad_linalg_test)

# --- mesh layer (node/edge/control-volume graph, tensor-grid producer) -------------
add_library(nitcad_mesh STATIC
  src/mesh/mesh.cpp
  src/mesh/tensor_grid.cpp)
add_library(NiTCAD::mesh ALIAS nitcad_mesh)
target_include_directories(nitcad_mesh PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(nitcad_mesh PUBLIC NiTCAD::base)

add_executable(nitcad_mesh_test
  tests/mesh/mesh_test.cpp
  tests/mesh/tensor_grid_test.cpp)
target_link_libraries(nitcad_mesh_test PRIVATE NiTCAD::mesh Catch2::Catch2WithMain)
add_test(NAME mesh COMMAND nitcad_mesh_test)

# --- physics layer (silicon parameters, Caughey-Thomas, SRH, Boltzmann statistics) --
add_library(nitcad_physics STATIC
  src/physics/semiconductor.cpp
  src/physics/mobility.cpp
  src/physics/recombination.cpp
  src/physics/statistics.cpp)
add_library(NiTCAD::physics ALIAS nitcad_physics)
target_include_directories(nitcad_physics PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(nitcad_physics PUBLIC NiTCAD::base)

add_executable(nitcad_physics_test
  tests/physics/semiconductor_test.cpp
  tests/physics/mobility_test.cpp
  tests/physics/recombination_test.cpp
  tests/physics/statistics_test.cpp)
target_link_libraries(nitcad_physics_test PRIVATE NiTCAD::physics Catch2::Catch2WithMain)
add_test(NAME physics COMMAND nitcad_physics_test)

# --- device layer (regions, doping, contacts; validated problem description) -------
add_library(nitcad_device STATIC src/device/device.cpp)
add_library(NiTCAD::device ALIAS nitcad_device)
target_include_directories(nitcad_device PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(nitcad_device PUBLIC NiTCAD::base NiTCAD::mesh NiTCAD::physics)

add_executable(nitcad_device_test
  tests/device/device_test.cpp
  tests/device/topology_test.cpp)
target_link_libraries(nitcad_device_test PRIVATE NiTCAD::device Catch2::Catch2WithMain)
add_test(NAME device COMMAND nitcad_device_test)

# --- assemble layer (scaling, Bernoulli/SG flux, ohmic values, equilibrium Poisson) --
add_library(nitcad_assemble STATIC
  src/assemble/scaling.cpp
  src/assemble/equilibrium_poisson.cpp)
add_library(NiTCAD::assemble ALIAS nitcad_assemble)
target_include_directories(nitcad_assemble PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(nitcad_assemble PUBLIC
  NiTCAD::base NiTCAD::linalg NiTCAD::mesh NiTCAD::physics NiTCAD::device)

add_executable(nitcad_assemble_test
  tests/assemble/flux_test.cpp
  tests/assemble/poisson_test.cpp)
target_link_libraries(nitcad_assemble_test PRIVATE NiTCAD::assemble Catch2::Catch2WithMain)
add_test(NAME assemble COMMAND nitcad_assemble_test)

# --- solve layer (Newton, equilibrium solve) ---------------------------------------
add_library(nitcad_solve STATIC src/solve/equilibrium.cpp)
add_library(NiTCAD::solve ALIAS nitcad_solve)
target_include_directories(nitcad_solve PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(nitcad_solve PUBLIC
  NiTCAD::base NiTCAD::linalg NiTCAD::assemble NiTCAD::device)

add_executable(nitcad_solve_test
  tests/solve/newton_test.cpp
  tests/solve/equilibrium_test.cpp)
target_link_libraries(nitcad_solve_test PRIVATE NiTCAD::solve Catch2::Catch2WithMain)
add_test(NAME solve COMMAND nitcad_solve_test)
