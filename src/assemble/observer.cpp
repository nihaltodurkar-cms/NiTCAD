#include "NiTCAD/assemble/observer.hpp"

namespace NiTCAD::assemble {

namespace detail {
Observer observer = nullptr;
}  // namespace detail

void set_observer(Observer observer) noexcept { detail::observer = observer; }

}  // namespace NiTCAD::assemble
