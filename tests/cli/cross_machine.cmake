# Runs a deck on this machine and compares its sweep with a committed reference made on the
# development machine (DECISIONS Unit 60, XC4-XC5; ARCHITECTURE 6.8: agreement by tolerances).
# -DNITCAD=<nitcad.exe> -DDECK=<deck> -DREFERENCE=<.nitsweep> -DWORK=<directory> -DEXPECT=<status>;
# optional -DREPLACE_FROM=<text> -DREPLACE_TO=<text>: the deck is edited first (a changed input,
# for the "differ" status, XG3).
file(REMOVE_RECURSE "${WORK}")
file(MAKE_DIRECTORY "${WORK}")
file(READ "${DECK}" text)
if(DEFINED REPLACE_FROM)
  string(FIND "${text}" "${REPLACE_FROM}" at)
  if(at EQUAL -1)
    message(FATAL_ERROR "the deck has no '${REPLACE_FROM}' to replace")
  endif()
  string(REPLACE "${REPLACE_FROM}" "${REPLACE_TO}" text "${text}")
endif()
get_filename_component(name "${DECK}" NAME)
file(WRITE "${WORK}/${name}" "${text}")
if(NOT text MATCHES "\n[ \t]*sweep[ \t]*=[ \t]*([^\r\n]+)")
  message(FATAL_ERROR "the deck writes no sweep")
endif()
string(STRIP "${CMAKE_MATCH_1}" sweep)

execute_process(COMMAND "${NITCAD}" run "${WORK}/${name}" RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "nitcad run: exit status ${rc}\n${out}${err}")
endif()
execute_process(COMMAND "${NITCAD}" compare "${REFERENCE}" "${WORK}/${sweep}"
                RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)
message("${out}${err}")
if(NOT rc EQUAL EXPECT)
  message(FATAL_ERROR "nitcad compare: exit status ${rc}, expected ${EXPECT}")
endif()
