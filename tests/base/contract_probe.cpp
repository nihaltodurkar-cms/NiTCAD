// Probe for NITCAD_EXPECTS. A violation terminates the process, so contract_test.cpp runs this
// executable as a child process: `contract_probe pass` must exit 0; `contract_probe violate` must
// terminate abnormally and print the violation text to stderr.
#include <windows.h>

#include <cstdio>
#include <cstring>

#include "NiTCAD/base/contract.hpp"

int main(int argc, char** argv) {
    // Keep the intentional crash from raising a Windows error-reporting dialog or delay.
    SetErrorMode(SEM_FAILCRITICALERRORS | SEM_NOGPFAULTERRORBOX);
    const bool violate = argc > 1 && std::strcmp(argv[1], "violate") == 0;
    // argc is not a compile-time constant, so neither branch is optimised away.
    NITCAD_EXPECTS(violate ? argc == 99 : argc >= 1);
    std::puts("contract satisfied");
    return 0;
}
