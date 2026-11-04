// VR2 (Unit 55, TVR-1; TVD-19): the reference library's public headers compile with only the include
// paths a nitcad_verify_reference consumer gets, which carry no Eigen. The include graph itself is checked
// by tests/verify/boundary.cmake before the library builds.
#include <catch2/catch_test_macros.hpp>

#include "NiTCAD/verify/reference/double_double.hpp"
#include "NiTCAD/verify/reference/gauss_kronrod.hpp"
#include "NiTCAD/verify/reference/tanh_sinh.hpp"

#if __has_include(<Eigen/Core>)
#error "Eigen's include path reaches nitcad_verify_reference consumers"
#endif

#ifdef EIGEN_WORLD_VERSION
#error "a verify/reference header includes Eigen"
#endif

TEST_CASE("verify/reference: public headers compile without Eigen") { SUCCEED("checked at compile time"); }
