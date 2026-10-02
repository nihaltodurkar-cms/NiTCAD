// Caughey-Thomas mobility gates (ARCHITECTURE.md section 11, Unit 5). Reference values were
// computed independently at 40 digits from the legacy formula and silicon parameters.
#include <catch2/catch_test_macros.hpp>

#include <cmath>

#include "NiTCAD/physics/mobility.hpp"

using namespace NiTCAD::physics;

namespace {

bool close(double a, double b, double rel) { return std::abs(a - b) <= rel * std::abs(b); }

struct Point {
    double N;
    double mu_n;
    double mu_p;
};

}  // namespace

TEST_CASE("mobility: legacy published-value gate for silicon at 300 K") {
    // Legacy test_caughey_thomas_matches_published_silicon_values.
    const Semiconductor si = silicon();
    REQUIRE(close(caughey_thomas_mobility(si, Carrier::electron, 0.0, 300.0), 1360.0, 1e-9));
    const double mu_1e18 = caughey_thomas_mobility(si, Carrier::electron, 1e18, 300.0);
    REQUIRE(close(mu_1e18, 300.0, 0.25));  // published ~270-300 cm^2/(V s)
}

TEST_CASE("mobility: silicon values across doping at 300 K") {
    const Semiconductor si = silicon();
    const Point table[] = {
        {0.0, 1360.0, 495.0},
        {1e15, 1345.0622836780626175, 476.59839908903143508},
        {1e16, 1247.9879239377467887, 406.43290552413809067},
        {1e17, 801.32674257569405969, 232.48156078201518299},
        {1e18, 263.30645839319100749, 96.450701225872366047},
        {1e19, 115.90787581787263789, 57.010279846043862874},
        {1e20, 94.990764017766295192, 49.346252246770342404},
    };
    for (const Point& pt : table) {
        CAPTURE(pt.N);
        REQUIRE(close(caughey_thomas_mobility(si, Carrier::electron, pt.N, 300.0), pt.mu_n, 1e-13));
        REQUIRE(close(caughey_thomas_mobility(si, Carrier::hole, pt.N, 300.0), pt.mu_p, 1e-13));
    }
}

TEST_CASE("mobility: temperature scales mu_max only") {
    const Semiconductor si = silicon();
    REQUIRE(close(caughey_thomas_mobility(si, Carrier::electron, 1e17, 400.0),
                  429.72201752612565687, 1e-13));
    REQUIRE(close(caughey_thomas_mobility(si, Carrier::hole, 1e17, 400.0),
                  135.65422348940712254, 1e-13));
    // At very high doping the mobility tends to mu_min, which does not depend on T.
    for (const double T : {250.0, 300.0, 400.0}) {
        REQUIRE(close(caughey_thomas_mobility(si, Carrier::electron, 1e30, T), 92.0, 1e-6));
    }
}

TEST_CASE("mobility: decreases monotonically with total impurity") {
    const Semiconductor si = silicon();
    for (const Carrier c : {Carrier::electron, Carrier::hole}) {
        double previous = caughey_thomas_mobility(si, c, 0.0, 300.0);
        for (double N = 1e12; N <= 1e22; N *= 1.5) {
            const double mu = caughey_thomas_mobility(si, c, N, 300.0);
            REQUIRE(mu < previous);
            previous = mu;
        }
    }
}
