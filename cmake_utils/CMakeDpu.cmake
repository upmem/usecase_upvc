if (NOT DEFINED DPU_DIRECTORY)
    message(FATAL_ERROR "You must set DPU_DIRECTORY variable")
endif()

SET(UPMEM_HOME $ENV{UPMEM_HOME})
SET(DPU_KCONFIG ${UPMEM_HOME}/usr/bin/dpukconfig)
SET(DPU_CC ${UPMEM_HOME}/usr/bin/dpucc)

macro(add_dpu_executable target configuration sources)
    dpu_generate_configuration(${target}_config.json ${configuration})
    dpu_build(${target} ${target}.bin ${target}_config.json "${sources}")
endmacro()

macro(dpu_generate_configuration output script)
    add_custom_command(
            OUTPUT ${CMAKE_BINARY_DIR}/${DPU_DIRECTORY}/${script}
            COMMAND ${CMAKE_COMMAND} -Dsource=${CMAKE_SOURCE_DIR}/${DPU_DIRECTORY} -Dscript=${script} -Ddestination=${CMAKE_BINARY_DIR}/${DPU_DIRECTORY} -Doutput=${output} -P ${CMAKE_SOURCE_DIR}/${CMAKE_UTILS_DIRECTORY}/CMakeDpuPrepareConfigurationScript.cmake
            DEPENDS ${script}
            VERBATIM
    )
    add_custom_command(
            OUTPUT ${output}
            COMMAND ${DPU_KCONFIG} -s ${CMAKE_BINARY_DIR}/${DPU_DIRECTORY}/${script}
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/${DPU_DIRECTORY}
            DEPENDS ${CMAKE_BINARY_DIR}/${DPU_DIRECTORY}/${script}
            VERBATIM
    )
endmacro()

macro(dpu_build target output configuration dpucc_opt sources)
    INCLUDE_DIRECTORIES(SYSTEM ${UPMEM_HOME}/usr/share/upmem/include/stdlib ${UPMEM_HOME}/usr/share/upmem/include/syslib)
    add_custom_target(${target} DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${output} SOURCES ${sources})
    get_property(incs TARGET ${target} PROPERTY INCLUDE_DIRECTORIES)
    foreach(inc ${incs})
        if (NOT ${inc} MATCHES "^${UPMEM_HOME}.*")
            set(__${target}_include_directories ${__${target}_include_directories} -I${inc})
        endif ()
    endforeach()

    add_executable(__dummy_exec_${target} EXCLUDE_FROM_ALL ${incs} ${sources})

    add_custom_command(
            OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${output}
            COMMAND ${DPU_CC} -b -C ${CMAKE_CURRENT_BINARY_DIR}/${configuration} ${dpucc_opt} -o ${CMAKE_CURRENT_BINARY_DIR}/${output} ${__${target}_include_directories} ${sources}
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/${DPU_DIRECTORY}
            DEPENDS ${configuration} ${sources}
            VERBATIM
    )
endmacro()

macro(add_dpu_dependency target dependency)
    add_dependencies(${target} ${dependency})
#    add_definitions(-DDPU_BINARY=${CMAKE_BINARY_DIR}/${DPU_DIRECTORY}/${dependency}.bin)
endmacro()