// The bundled examples (Unit 26, VQ11, WG8): exactly the four decks, titled and described, each running
// to completion; a user's copy kept.
#include <catch2/catch_test_macros.hpp>

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "NiTCAD/app/examples.hpp"
#include "NiTCAD/app/run_controller.hpp"
#include "fixtures.hpp"

using namespace NiTCAD;
namespace fs = std::filesystem;

namespace {

const fs::path source_examples{NITCAD_EXAMPLES_DIR};

}  // namespace

TEST_CASE("app examples: the four decks, titled and described (WG8)", "[app]") {
    const auto examples = app::list_examples(source_examples);
    REQUIRE(examples.size() == 4);
    const std::vector<std::pair<std::string, std::string>> expected{{"moscap_cv.deck", "MOS capacitor"},
                                                                    {"mosfet_2d.deck", "2D MOSFET"},
                                                                    {"mosfet_3d.deck", "3D MOSFET"},
                                                                    {"pn_diode_1d.deck", "1D pn diode"}};
    for (std::size_t k = 0; k < 4; ++k) {
        CAPTURE(k);
        CHECK(examples[k].path.filename() == expected[k].first);
        CHECK(examples[k].title == expected[k].second);
        CHECK_FALSE(examples[k].description.empty());
    }
    CHECK(examples[3].description == "A silicon pn junction swept forward: the band diagram and the I-V curve.");
    // A directory without decks has none; a deck without a comment is titled by its name.
    const fs::path dir = app_test::scratch("examples_list");
    CHECK(app::list_examples(dir).empty());
    CHECK(app::list_examples(dir / "missing").empty());
    app_test::write(dir, "plain.deck", "[mesh]\n");
    const auto plain = app::list_examples(dir);
    REQUIRE(plain.size() == 1);
    CHECK(plain[0].title == "plain");
    CHECK(plain[0].description.empty());
    // Beside the executable when there, else the fallback.
    CHECK(app::examples_directory(dir / "nowhere", source_examples) == source_examples);
    fs::create_directories(dir / "bin" / "examples");
    app_test::write(dir / "bin" / "examples", "a.deck", "# A\n");
    CHECK(app::examples_directory(dir / "bin", source_examples) == dir / "bin" / "examples");
}

TEST_CASE("app examples: a copy to run, kept once made", "[app]") {
    const fs::path to = app_test::scratch("examples_copy") / "mine";
    const auto examples = app::list_examples(source_examples);
    REQUIRE_FALSE(examples.empty());
    auto copy = app::copy_example(examples[0], to);
    REQUIRE(copy);
    CHECK(*copy == to / examples[0].path.filename());
    CHECK(fs::exists(*copy));
    // The user's edits survive opening the example again.
    std::ofstream(*copy, std::ios::binary | std::ios::app) << "# edited\n";
    const auto size = fs::file_size(*copy);
    auto again = app::copy_example(examples[0], to);
    REQUIRE(again);
    CHECK(fs::file_size(*again) == size);
}

TEST_CASE("app examples: every example deck runs to completion (WG8)", "[app]") {
    const fs::path dir = app_test::scratch("examples_run");
    for (const app::Example& e : app::list_examples(source_examples)) {
        CAPTURE(e.path.filename().string());
        auto copy = app::copy_example(e, dir);
        REQUIRE(copy);
        std::ifstream f(*copy, std::ios::binary);
        const std::string text{std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>()};
        app::RunController rc;
        REQUIRE(rc.start_run(text, *copy));
        rc.wait();
        rc.process();
        INFO((rc.error() ? rc.error()->message : std::string()));
        REQUIRE(rc.state() == app::RunState::done);
        CHECK_FALSE(rc.stopped_early());
        const auto& r = *rc.report();
        const std::size_t points = r.sweep ? r.sweep->points.size() : r.small_signal ? r.small_signal->points.size() : 0;
        CHECK(points == rc.plan()->points.size());
        CHECK(points >= 5);
        CHECK(r.written.size() >= 3);  // the result, the device, the table
    }
}
