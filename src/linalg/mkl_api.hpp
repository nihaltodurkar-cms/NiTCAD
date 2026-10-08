// The MKL entry points the PARDISO backend calls, resolved from the runtime that
// load_mkl_runtime loaded (mkl_runtime.hpp). Private to linalg. The signatures are MKL's (oneMKL
// developer reference) with MKL_INT = int, the LP64 interface layer the loader selects.
#pragma once

namespace NiTCAD::linalg::detail {

using MklInt = int;

struct MklApi {
    void (*pardisoinit)(void* pt, const MklInt* mtype, MklInt* iparm);
    void (*pardiso)(void* pt, const MklInt* maxfct, const MklInt* mnum, const MklInt* mtype,
                    const MklInt* phase, const MklInt* n, const void* a, const MklInt* ia,
                    const MklInt* ja, MklInt* perm, const MklInt* nrhs, MklInt* iparm,
                    const MklInt* msglvl, void* b, void* x, MklInt* error);
    int (*set_num_threads_local)(int threads);  // MKL_Set_Num_Threads_Local: thread-local
};

// The entry points of the loaded runtime; nullptr until load_mkl_runtime has succeeded.
[[nodiscard]] const MklApi* mkl_api() noexcept;

}  // namespace NiTCAD::linalg::detail
