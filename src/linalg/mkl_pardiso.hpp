// Intel MKL PARDISO backend (Unit 18). Private to linalg; MKL is reached only through the entry
// points of the runtime load_mkl_runtime loaded (mkl_api.hpp), never through an MKL header.
#pragma once

#include <array>
#include <complex>
#include <cstddef>
#include <expected>
#include <span>
#include <type_traits>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/linalg/sparse_matrix.hpp"
#include "backend.hpp"
#include "mkl_api.hpp"

namespace NiTCAD::linalg::detail {

// Real unsymmetric (matrix type 11) or complex unsymmetric (13) PARDISO, on the CSR pattern as
// given (zero-based, the layout PARDISO reads natively). The analysis (phase 11: weighted matching
// and scaling, which need values, then METIS nested dissection and the symbolic factorization) runs
// at the first factorization after analyze(); every later factorization of the pattern is phase 22
// only, every solve phase 33. Requires the MKL runtime to be loaded.
template <MatrixScalar Scalar>
class PardisoBackend final : public Backend<Scalar> {
public:
    explicit PardisoBackend(int threads);
    ~PardisoBackend() override;
    PardisoBackend(const PardisoBackend&) = delete;
    PardisoBackend& operator=(const PardisoBackend&) = delete;

    [[nodiscard]] std::expected<void, base::Error> analyze(
        const BasicSparseMatrix<Scalar>& a) override;
    [[nodiscard]] std::expected<Factorization, base::Error> factorize(
        std::span<const Scalar> csr_values) override;
    // A failed phase 33 leaves x NaN, which BasicLinearSolver reports as inaccurate_solve.
    void solve(std::span<const Scalar> b, std::span<Scalar> x) override;

private:
    // One PARDISO call with this solver's thread count set thread-locally around it; counts the
    // analysis, factorization and solve phases it contains. Returns PARDISO's error code.
    MklInt call(MklInt phase, Scalar* b, Scalar* x);
    void release() noexcept;  // phase -1, if anything is held

    static constexpr MklInt matrix_type = std::is_same_v<Scalar, double> ? 11 : 13;

    const MklApi* api_;
    int threads_;
    std::array<void*, 64> pt_{};  // PARDISO's handle: must not move while it holds a factorization
    std::array<MklInt, 64> iparm_{};
    MklInt n_ = 0;
    std::vector<MklInt> row_offsets_;
    std::vector<MklInt> columns_;
    std::vector<Scalar> values_;
    std::vector<Scalar> rhs_;  // PARDISO's b is not const
    bool held_ = false;      // pt_ holds an analysis
    bool analyzed_ = false;  // phase 11 done for the current pattern
};

}  // namespace NiTCAD::linalg::detail
