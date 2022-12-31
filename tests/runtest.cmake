message("Running: ${TESTCMD} with input ${TESTINP}")
execute_process(COMMAND ${TESTCMD} INPUT_FILE ${TESTINP} OUTPUT_FILE ${TESTOUT} RESULT_VARIABLE RET)

file(READ ${TESTOUT} TEST_OUTPUT)
message("Test output: ${TEST_OUTPUT}")

if(NOT (RET EQUAL 0))
  message(FATAL_ERROR "Test program ${TESTCMD} returned ${RET}")
endif()

