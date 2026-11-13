// `nitcad learn` (Unit 17b; DECISIONS.md Unit 17 ID8, ID10, Unit 17b LB6, LB7):
//
//   nitcad learn manifest <manifest> <history.nithist>...
//       loads each adaptation history (refusing a corrupt one or another schema), digests it and
//       writes the manifest, sorted by digest (io/learn_files.hpp).
//
// Exit status: 0 done; 1 a file failed (the message on stderr); 64 usage.
#pragma once

namespace NiTCAD::cli {

// argv[1] is "learn".
int learn(int argc, char** argv);

}  // namespace NiTCAD::cli
