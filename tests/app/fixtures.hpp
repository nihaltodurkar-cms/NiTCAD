// Shared fixtures of the app tests (Unit 26): decks written to scratch directories, and running the
// nitcad executable.
#pragma once

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>

#include <catch2/catch_test_macros.hpp>

namespace app_test {

inline std::filesystem::path scratch(const std::string& name) {
    const std::filesystem::path dir = std::filesystem::temp_directory_path() / "nitcad_app_test" / name;
    std::filesystem::remove_all(dir);
    std::filesystem::create_directories(dir);
    return dir;
}

inline std::filesystem::path write(const std::filesystem::path& dir, const std::string& name, const std::string& text) {
    const std::filesystem::path p = dir / name;
    std::ofstream(p, std::ios::binary) << text;
    return p;
}

// A 1D pn diode (the deck tests' diode) with its outputs.
inline const std::string diode_device = R"(# A 1D pn diode
[mesh]
type = tensor
x = linspace(0, 2um, 201)

[region silicon]
material = silicon

[doping]
profile = acceptor constant 1e17 x=-1:1um
profile = donor constant 1e17 x=1um:3um

[contact anode]
kind = ohmic
patch = x_min

[contact cathode]
kind = ohmic
patch = x_max
)";

inline std::string diode_sweep(const std::string& series = "0 : 0.4 : 0.1") {
    return diode_device + "\n[sweep]\ncontact anode = " + series +
           "\n\n[output]\nsweep = diode.nitsweep\ndevice = diode.nitdev\ntable = diode.csv\n";
}

// A small 2D n-channel MOSFET on a tensor mesh: n+ source and drain pockets under a 50 nm oxide, a
// metal gate over the channel, the body below; the gate swept at a drain bias of 0.5 V.
inline const std::string mosfet2d = R"(# A 2D n-channel MOSFET
[mesh]
type = tensor
x = linspace(0, 2um, 21)
y = linspace(0, 1um, 21)

[region silicon]
material = silicon

[region oxide]
material = sio2
box = 0 2um  0.95um 1um

[doping]
profile = acceptor constant 1e16
profile = donor gaussian 2e19 axis=y center=0.95um width=0.15um x=0:0.5um
profile = donor gaussian 2e19 axis=y center=0.95um width=0.15um x=1.5um:2um

[contact source]
kind = ohmic
patch = x_min
box = 0 0  0.6um 0.9um

[contact drain]
kind = ohmic
patch = x_max
box = 2um 2um  0.6um 0.9um

[contact gate]
kind = electrode
patch = y_max
box = 0.6um 1.4um  1um 1um
electrode = metal
work_function = 4.1

[contact body]
kind = ohmic
patch = y_min

[sweep]
contact gate = 0 : 1.5 : 0.25
contact drain = 0.5

[output]
sweep = mosfet.nitsweep
device = mosfet.nitdev
)";

// Generated 2D and 3D resistors at equilibrium, and the diode's transient.
inline const std::string generate2d = R"(# A generated 2D resistor
[mesh]
type = generate2d
patch left = 0 0, 0 1um
patch right = 2um 0, 2um 1um
max_size = 0.25um

[region silicon]
material = silicon
polygon = 0 0, 2um 0, 2um 1um, 0 1um

[doping]
profile = donor constant 1e16

[contact left]
kind = ohmic
patch = left

[contact right]
kind = ohmic
patch = right

[equilibrium]
)";

inline const std::string generate3d = R"(# A generated 3D resistor
[mesh]
type = generate3d
cuboid = 0 1um 0 1um 0 1um
max_size = 0.5um

[region silicon]
material = silicon

[doping]
profile = donor constant 1e16

[contact left]
kind = ohmic
patch = x_min

[contact right]
kind = ohmic
patch = x_max

[equilibrium]
)";

inline std::string transient_diode() {
    return diode_device + "\n[transient]\ncontact anode = step(0, 0.3, 1e-9)\nt_end = 4e-9\n";
}

struct Cli {
    int status;
    std::string output;
};

// Runs the nitcad executable with these arguments, stdout and stderr together.
inline Cli nitcad(const std::string& arguments) {
    const std::string command = "\"\"" NITCAD_CLI_PATH "\" " + arguments + " 2>&1\"";
    FILE* pipe = _popen(command.c_str(), "r");
    REQUIRE(pipe != nullptr);
    std::string out;
    char buffer[512];
    while (std::fgets(buffer, sizeof buffer, pipe)) out += buffer;
    return {_pclose(pipe), out};
}

inline std::string quoted(const std::filesystem::path& p) { return "\"" + p.string() + "\""; }

}  // namespace app_test
