# Writes OUTPUT, the build facts of a run's provenance (Unit 42; DECISIONS.md PV3), at every build:
# the version, build type and compiler given by the configure step, and the source revision and
# dirty state read from git now. The file is replaced only when its text changes, so an unchanged
# revision rebuilds nothing. Outside a git checkout, or without git, the revision is "unknown".
#
#   cmake -DSOURCE=<dir> -DOUTPUT=<file> -DGIT=<git or empty> -DVERSION=<v> -DBUILD_TYPE=<t>
#         -DCOMPILER=<c> -P build_info.cmake
set(revision "unknown")
set(dirty "false")
if(GIT)
  execute_process(COMMAND "${GIT}" -C "${SOURCE}" rev-parse HEAD
                  OUTPUT_VARIABLE head RESULT_VARIABLE failed
                  OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
  if(NOT failed)
    set(revision "${head}")
    # Untracked files count: a new source file not yet committed changes what is built.
    execute_process(COMMAND "${GIT}" -C "${SOURCE}" status --porcelain
                    OUTPUT_VARIABLE changes RESULT_VARIABLE failed ERROR_QUIET)
    if(failed OR NOT changes STREQUAL "")
      set(dirty "true")
    endif()
  endif()
endif()

set(text "// Generated at build time by cmake/build_info.cmake (Unit 42). Do not edit.
#pragma once

namespace NiTCAD::build_info {

inline constexpr const char* version = \"${VERSION}\";
inline constexpr const char* revision = \"${revision}\";
inline constexpr bool dirty = ${dirty};
inline constexpr const char* build_type = \"${BUILD_TYPE}\";
inline constexpr const char* compiler = \"${COMPILER}\";

}  // namespace NiTCAD::build_info
")
file(WRITE "${OUTPUT}.tmp" "${text}")
file(COPY_FILE "${OUTPUT}.tmp" "${OUTPUT}" ONLY_IF_DIFFERENT)
file(REMOVE "${OUTPUT}.tmp")
