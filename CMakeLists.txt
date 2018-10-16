project(upvc)
set( CMAKE_C_FLAGS "-O3")
add_executable(upvc code.c compare.c dispatch.c getgenome.c getread.c index.c processread.c vcf.c upvc_dpu.c upvc_host.c vartree.c)
target_link_libraries(upvc pthread)
