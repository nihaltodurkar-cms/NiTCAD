// The wall time of a public run, for its provenance (Unit 42; solve/provenance.hpp): from the call
// to its result, on the monotonic clock.
#pragma once

#include <chrono>

namespace NiTCAD::solve::detail {

class WallClock {
public:
    [[nodiscard]] double seconds() const {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - start_).count();
    }

private:
    std::chrono::steady_clock::time_point start_ = std::chrono::steady_clock::now();
};

}  // namespace NiTCAD::solve::detail
