// Unit 68 slice 1 (MAC9; DECISIONS.md amendment 5, S10): the SiGe-base HBT library deck runs its Gummel
// plot, and its collector current at V_BE = 0.7 V exceeds the all-silicon base's by the textbook factor
// (mu_n ni_eff^2 of the base, which carries exp(dEg / kT) with the mobility, Nc Nv and BGN changes) within 20 %.
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>

#include "NiTCAD/deck/run.hpp"
#include "NiTCAD/physics/bandgap_narrowing.hpp"
#include "NiTCAD/physics/mobility.hpp"
#include "NiTCAD/physics/sige.hpp"
#include "fixtures.hpp"

using namespace NiTCAD;

namespace {

std::string hbt_deck() {
    std::ifstream in(std::filesystem::path(NITCAD_EXAMPLES_DIR) / "devices" / "sige_hbt_2d.deck");
    std::stringstream s;
    s << in.rdbuf();
    return s.str();
}

// The collector current at the sweep's last point (V_BE = V_CB + 0.7 V = 0.7 V).
double collector_current(const std::string& text, const std::string& name) {
    const deck::Plan p = deck_test::must(text, deck_test::scratch(name));
    const auto r = deck::run(p);
    if (!r) FAIL(r.error().message);
    REQUIRE(r->sweep.has_value());
    REQUIRE_FALSE(r->sweep->stopped.has_value());
    const auto& last = r->sweep->points.back();
    REQUIRE(std::abs(last.bias_V[1] - 0.7) < 1e-12);
    return last.terminal_current[0];
}

double base_factor(const physics::SemiconductorParameters& p) {
    const auto m = *physics::Semiconductor::create(p);
    const double ni = physics::effective_intrinsic_density(m, 2e18, 300.0);
    return physics::caughey_thomas_mobility(m, physics::Carrier::electron, 2e18, 300.0) * ni * ni;
}

}  // namespace

TEST_CASE("MAC9: the SiGe HBT's collector current against the all-silicon base", "[deck][unit68][sige]") {
    const std::string sige = hbt_deck();
    REQUIRE(sige.find("xge = 0.2  #") != std::string::npos);
    std::string silicon = sige;
    silicon.replace(silicon.find("xge = 0.2  #"), 11, "xge = 0.0  #");
    const double ic_sige = collector_current(sige, "mac9_sige");
    const double ic_si = collector_current(silicon, "mac9_si");
    const double ratio = ic_sige / ic_si;
    const double textbook = base_factor(*physics::sige_parameters(0.2, true)) / base_factor(physics::silicon_parameters);
    const double dEg = physics::sige_bands(0.2, 0.0)->gap_change;
    std::printf("MAC9 Ic(SiGe) %.4e  Ic(Si) %.4e A/cm  ratio %.1f  textbook %.1f (exp(-dEg/kT) alone %.1f)  off by %+.1f %%\n",
                ic_sige, ic_si, ratio, textbook, std::exp(-dEg / 0.025852), 100.0 * (ratio / textbook - 1.0));
    CHECK(std::abs(ratio / textbook - 1.0) <= 0.20);
}
