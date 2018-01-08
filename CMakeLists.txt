project(genomee)

add_executable(genomee code.c compare.c dispatch.c getgenome.c getread.c index.c processread.c vcf.c upvc_dpu.c upvc_host.c)
target_link_libraries(genomee pthread)
