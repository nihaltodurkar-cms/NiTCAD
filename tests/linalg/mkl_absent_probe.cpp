// A process that never loads MKL (Unit 18): mkl_pardiso solvers must refuse to start, and a failed
// load must leave nothing loaded. Its own executable, because the PARDISO tests load the runtime
// and Catch2 runs them in random order. Exit 0 when every check holds, else the failed check's
// number.
#include "NiTCAD/linalg/linear_solver.hpp"
#include "NiTCAD/linalg/mkl_runtime.hpp"

using namespace NiTCAD::linalg;
using NiTCAD::base::ErrorCode;

int main() {
    const SolverConfig pardiso{.backend = SolverBackend::mkl_pardiso};
    if (loaded_mkl_runtime()) return 1;
    const auto real = LinearSolver::create(pardiso);
    if (real || real.error().code != ErrorCode::invalid_input) return 2;
    const auto complex = ComplexLinearSolver::create(pardiso);
    if (complex || complex.error().code != ErrorCode::invalid_input) return 3;
    const auto bad = load_mkl_runtime("no/such/mkl_rt.dll");
    if (bad || bad.error().code != ErrorCode::invalid_input) return 4;
    if (loaded_mkl_runtime()) return 5;
    if (LinearSolver::create(pardiso)) return 6;
    if (!LinearSolver::create({})) return 7;  // Eigen needs nothing
    return 0;
}
