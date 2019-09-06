cmake_minimum_required(VERSION 2.8.12)

add_subdirectory(host)
add_subdirectory(tests)

add_custom_target(plotdata ${CMAKE_CURRENT_BINARY_DIR}/tests/upvclog2graph ${CMAKE_CURRENT_SOURCE_DIR}/tests/chr22_integration/upvc_log.txt
	COMMAND echo "cycles" > plotdata.csv
	COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests/moy_time.bash ${CMAKE_CURRENT_SOURCE_DIR}/tests/chr22_integration/upvc_log.txt.csv | sed -e "'s/^.* \\([^ ]*\\)$$/\\1/'" | head -n 1 >> plotdata.csv
	COMMAND cat plotdata.csv)
