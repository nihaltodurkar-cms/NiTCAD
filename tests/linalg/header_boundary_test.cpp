// Backend neutrality check (ARCHITECTURE.md 6.10, Q5). This translation unit includes every
// public linalg header and is compiled with the include paths a linalg consumer gets. It fails
// to build if Eigen's include path reaches consumers, or if a public header includes Eigen by
// any path.
#include <catch2/catch_test_macros.hpp>

#include "NiTCAD/linalg/linear_solver.hpp"
#include "NiTCAD/linalg/sparse_matrix.hpp"

#if __has_include(<Eigen/Core>)
#error "Eigen's include path reaches linalg consumers: NiTCAD::linalg must link Eigen privately"
#endif

#ifdef EIGEN_WORLD_VERSION
#error "a public linalg header includes Eigen"
#endif

TEST_CASE("linalg: public headers compile without Eigen") {
    SUCCEED("checked at compile time");
}
