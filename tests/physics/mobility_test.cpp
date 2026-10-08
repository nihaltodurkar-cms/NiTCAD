// Caughey-Thomas mobility gates (ARCHITECTURE.md section 11, Unit 5). References computed here in
// double-double arithmetic (references.hpp) from the legacy formula and silicon parameters.
#include <catch2/catch_test_macros.hpp>

#include <cmath>

#include "NiTCAD/physics/mobility.hpp"
#include "references.hpp"

using namespace NiTCAD::physics;

namespace {

bool close(double a, double b, double rel) { return std::abs(a - b) <= rel * std::abs(b); }

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
    const auto& e = silicon_parameters.electron_mobility;
    const auto& h = silicon_parameters.hole_mobility;
    REQUIRE(caughey_thomas_mobility(si, Carrier::electron, 0.0, 300.0) == 1360.0);
    REQUIRE(caughey_thomas_mobility(si, Carrier::hole, 0.0, 300.0) == 495.0);
    for (const double N : {1e15, 1e16, 1e17, 1e18, 1e19, 1e20}) {
        CAPTURE(N);
        REQUIRE(close(caughey_thomas_mobility(si, Carrier::electron, N, 300.0),
                      reference::caughey_thomas(e, N, 300.0).value(), 1e-13));
        REQUIRE(close(caughey_thomas_mobility(si, Carrier::hole, N, 300.0),
                      reference::caughey_thomas(h, N, 300.0).value(), 1e-13));
    }
}

TEST_CASE("mobility: temperature scales mu_max only") {
    const Semiconductor si = silicon();
    REQUIRE(close(caughey_thomas_mobility(si, Carrier::electron, 1e17, 400.0),
                  reference::caughey_thomas(silicon_parameters.electron_mobility, 1e17, 400.0)
                      .value(),
                  1e-13));
    REQUIRE(close(caughey_thomas_mobility(si, Carrier::hole, 1e17, 400.0),
                  reference::caughey_thomas(silicon_parameters.hole_mobility, 1e17, 400.0).value(),
                  1e-13));
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
