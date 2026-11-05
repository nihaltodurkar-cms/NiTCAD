// VR2 (Unit 55, TVR-1; TVD-19): the consistency layer's public headers reach no physics, assemble, solve, linalg
// or device header, directly or transitively, and no Eigen (the test links the deck, whose include path carries
// Eigen, so only the macro is checked here; the reference layer checks the path). Each forbidden layer has one of its types declared
// here first as a variable of the same name: a header of that layer, included through any path, redeclares the
// name as a class and the translation unit does not compile. tests/verify/boundary.cmake checks the direct
// include graph before the library builds.
namespace NiTCAD::physics {
inline int Semiconductor = 0;
}
namespace NiTCAD::assemble {
inline int DriftDiffusion = 0;
}
namespace NiTCAD::solve {
inline int BiasOptions = 0;
}
namespace NiTCAD::linalg {
inline int SolverConfig = 0;
}
namespace NiTCAD::device {
inline int Device = 0;
}

#include "NiTCAD/verify/consistency/checks.hpp"
#include "NiTCAD/verify/consistency/device_data.hpp"


#ifdef EIGEN_WORLD_VERSION
#error "a verify/consistency header includes Eigen"
#endif

#include <catch2/catch_test_macros.hpp>

TEST_CASE("verify/consistency: public headers reach no forbidden layer", "[verify][boundary]") {
    CHECK(NiTCAD::physics::Semiconductor == 0);
    CHECK(NiTCAD::device::Device == 0);
}
