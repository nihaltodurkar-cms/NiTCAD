// `nitcad learn` (Unit 17b; learn.hpp).
#include "learn.hpp"

#include <cstdio>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

#include "NiTCAD/io/learn_files.hpp"

namespace NiTCAD::cli {
namespace {

int usage() {
    std::fputs("usage: nitcad learn manifest <manifest> <history.nithist>...\n", stderr);
    return 64;
}

int fail(const base::Error& e) {
    std::fprintf(stderr, "nitcad learn: %s\n", e.message.c_str());
    return 1;
}

int manifest(int argc, char** argv) {
    if (argc < 5) return usage();  // nitcad learn manifest <manifest> <history>...
    std::vector<std::filesystem::path> histories;
    for (int i = 4; i < argc; ++i) histories.emplace_back(argv[i]);
    auto m = io::make_manifest(histories);
    if (!m) return fail(m.error());
    if (auto w = io::save_manifest(argv[3], *m); !w) return fail(w.error());
    for (const results::ManifestEntry& e : m->entries) {
        std::printf("%016llx %s\n", static_cast<unsigned long long>(e.digest), e.path.c_str());
    }
    std::printf("manifest: %zu histories, schema %llu\n", m->entries.size(),
                static_cast<unsigned long long>(m->schema));
    return 0;
}

}  // namespace

int learn(int argc, char** argv) {
    if (argc < 3) return usage();
    const std::string_view sub = argv[2];
    if (sub == "manifest") return manifest(argc, argv);
    return usage();
}

}  // namespace NiTCAD::cli
