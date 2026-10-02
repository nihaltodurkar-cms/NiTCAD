#include <catch2/catch_test_macros.hpp>

#include <array>
#include <vector>

#include "NiTCAD/linalg/sparse_matrix.hpp"

using namespace NiTCAD::linalg;
using NiTCAD::base::ErrorCode;

namespace {
std::vector<Index> to_vector(std::span<const Index> s) { return {s.begin(), s.end()}; }
std::vector<double> to_vector(std::span<const double> s) { return {s.begin(), s.end()}; }
}  // namespace

TEST_CASE("sparse matrix: default is the empty 0 x 0 matrix") {
    const SparseMatrix m;
    REQUIRE(m.rows() == 0);
    REQUIRE(m.cols() == 0);
    REQUIRE(m.nonzeros() == 0);
    REQUIRE(to_vector(m.row_offsets()) == std::vector<Index>{0});
}

TEST_CASE("sparse matrix: triplets become sorted CSR with duplicates summed and zeros kept") {
    // [ 1  0  2 ]
    // [ 0  0  0 ]      row 1 is empty
    // [ 0  5  0 ]      plus an explicit zero at (2, 0)
    const std::array<Triplet, 6> entries{{
        {2, 1, 2.0},
        {0, 2, 2.0},
        {0, 0, 1.0},
        {2, 0, 0.0},
        {2, 1, 3.0},  // duplicate of (2, 1)
        {0, 0, 0.0},  // duplicate of (0, 0)
    }};
    const auto m = SparseMatrix::from_triplets(3, 3, entries);
    REQUIRE(m.has_value());
    REQUIRE(m->rows() == 3);
    REQUIRE(m->cols() == 3);
    REQUIRE(m->nonzeros() == 4);
    REQUIRE(to_vector(m->row_offsets()) == std::vector<Index>{0, 2, 2, 4});
    REQUIRE(to_vector(m->col_indices()) == std::vector<Index>{0, 2, 0, 1});
    REQUIRE(to_vector(m->values()) == std::vector<double>{1.0, 2.0, 0.0, 5.0});
}

TEST_CASE("sparse matrix: duplicates are summed in input order") {
    // 1e16 + 1 - 1e16 is 0 in this order and 1 in another, so the order is observable.
    const std::array<Triplet, 3> entries{{{0, 0, 1e16}, {0, 0, 1.0}, {0, 0, -1e16}}};
    const auto m = SparseMatrix::from_triplets(1, 1, entries);
    REQUIRE(m.has_value());
    REQUIRE(m->values()[0] == 0.0);
}

TEST_CASE("sparse matrix: rectangular matrices and empty entry lists are allowed") {
    const std::array<Triplet, 2> entries{{{1, 3, 4.0}, {0, 0, 1.0}}};
    const auto m = SparseMatrix::from_triplets(2, 4, entries);
    REQUIRE(m.has_value());
    REQUIRE(to_vector(m->row_offsets()) == std::vector<Index>{0, 1, 2});

    const auto empty = SparseMatrix::from_triplets(3, 3, {});
    REQUIRE(empty.has_value());
    REQUIRE(empty->nonzeros() == 0);
    REQUIRE(to_vector(empty->row_offsets()) == std::vector<Index>{0, 0, 0, 0});
}

TEST_CASE("sparse matrix: invalid dimensions and indices are errors, not clamped") {
    REQUIRE(SparseMatrix::from_triplets(-1, 2, {}).error().code == ErrorCode::invalid_input);
    REQUIRE(SparseMatrix::from_triplets(2, -1, {}).error().code == ErrorCode::invalid_input);

    const std::array<Triplet, 3> bad_row{{{0, 0, 1.0}, {1, 1, 1.0}, {2, 0, 1.0}}};
    const auto r = SparseMatrix::from_triplets(2, 2, bad_row);
    REQUIRE_FALSE(r.has_value());
    REQUIRE(r.error().code == ErrorCode::invalid_input);
    REQUIRE(r.error().context->index == 2);  // the offending entry

    const std::array<Triplet, 1> bad_col{{{0, -1, 1.0}}};
    const auto c = SparseMatrix::from_triplets(2, 2, bad_col);
    REQUIRE(c.error().code == ErrorCode::invalid_input);
    REQUIRE(c.error().context->index == 0);
}

TEST_CASE("sparse matrix: multiply computes y = A x") {
    // [ 2 -1  0 ]   [1]   [ 0]
    // [ 0  3  4 ] * [2] = [18]
    const std::array<Triplet, 4> entries{{{0, 0, 2.0}, {0, 1, -1.0}, {1, 1, 3.0}, {1, 2, 4.0}}};
    const auto m = SparseMatrix::from_triplets(2, 3, entries);
    REQUIRE(m.has_value());
    const std::array<double, 3> x{1.0, 2.0, 3.0};
    std::array<double, 2> y{-7.0, -7.0};
    m->multiply(x, y);
    REQUIRE(y == std::array<double, 2>{0.0, 18.0});
}

TEST_CASE("sparse matrix: values are rewritten in place and the pattern is unchanged") {
    const std::array<Triplet, 3> entries{{{0, 0, 1.0}, {0, 1, 2.0}, {1, 1, 3.0}}};
    auto m = SparseMatrix::from_triplets(2, 2, entries);
    REQUIRE(m.has_value());
    const SparseMatrix before = *m;
    for (double& v : m->values()) {
        v *= 10.0;
    }
    REQUIRE(to_vector(m->values()) == std::vector<double>{10.0, 20.0, 30.0});
    REQUIRE(m->has_same_pattern(before));
}

TEST_CASE("sparse matrix: pattern comparison sees dimensions, offsets and columns") {
    const std::array<Triplet, 2> diag{{{0, 0, 1.0}, {1, 1, 1.0}}};
    const std::array<Triplet, 2> moved{{{0, 1, 1.0}, {1, 1, 1.0}}};
    const std::array<Triplet, 3> more{{{0, 0, 1.0}, {0, 1, 1.0}, {1, 1, 1.0}}};
    const auto a = SparseMatrix::from_triplets(2, 2, diag);
    const auto same_values_changed = SparseMatrix::from_triplets(
        2, 2, std::array<Triplet, 2>{{{0, 0, 5.0}, {1, 1, -5.0}}});
    REQUIRE(a->has_same_pattern(*same_values_changed));
    REQUIRE_FALSE(a->has_same_pattern(*SparseMatrix::from_triplets(2, 2, moved)));
    REQUIRE_FALSE(a->has_same_pattern(*SparseMatrix::from_triplets(2, 2, more)));
    REQUIRE_FALSE(a->has_same_pattern(*SparseMatrix::from_triplets(3, 3, diag)));
}
