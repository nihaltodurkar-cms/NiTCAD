// The MKL runtime the build was configured with (CMake cache variable NITCAD_MKL_RUNTIME), loaded
// once for the PARDISO tests. REQUIRE_MKL() skips a test when none is configured and fails it when
// the configured one does not load.
#pragma once

#include <catch2/catch_test_macros.hpp>

#include <expected>
#include <filesystem>
#include <string>
#include <string_view>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/linalg/mkl_runtime.hpp"

namespace mkl_test {

inline constexpr std::string_view configured_path = NITCAD_MKL_RUNTIME_PATH;

inline const std::expected<NiTCAD::linalg::MklRuntime, NiTCAD::base::Error>& runtime() {
    static const auto loaded =
        NiTCAD::linalg::load_mkl_runtime(std::filesystem::path(std::string(configured_path)));
    return loaded;
}

}  // namespace mkl_test

#define REQUIRE_MKL()                                                                   \
    do {                                                                                \
        if (mkl_test::configured_path.empty()) {                                        \
            SKIP("no MKL runtime configured (NITCAD_MKL_RUNTIME)");                      \
        }                                                                               \
        INFO((mkl_test::runtime() ? std::string("MKL ") + mkl_test::runtime()->version  \
                                  : mkl_test::runtime().error().message));              \
        REQUIRE(mkl_test::runtime().has_value());                                       \
    } while (false)
