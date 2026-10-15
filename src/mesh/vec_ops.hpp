// Small vector algebra over any arithmetic type (double, the predicates' filtered doubles, exact
// integers), shared by the predicates and the dual measures. Internal to the mesh layer.
#pragma once

#include <array>

namespace NiTCAD::mesh::detail {

template <class T>
using Vec = std::array<T, 3>;

template <class T>
Vec<T> sub(const Vec<T>& a, const Vec<T>& b) {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}
template <class T>
T dot3(const Vec<T>& a, const Vec<T>& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
template <class T>
Vec<T> cross(const Vec<T>& a, const Vec<T>& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
template <class T>
Vec<T> scale(const T& s, const Vec<T>& a) {
    return {s * a[0], s * a[1], s * a[2]};
}
template <class T>
Vec<T> add(const Vec<T>& a, const Vec<T>& b) {
    return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

}  // namespace NiTCAD::mesh::detail
