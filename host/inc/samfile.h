#ifndef __SAMFILE_H__
#define __SAMFILE_H__
#include "genome.h"

#define SAM_VERSION "1.6"
#define PROGRAM_ID  "1"
#define PROGRAM_NAME "Mappim"
#define PROGRAM_DESCRIPTION "A mapper using processing in memory accelerations"
#define PROGRAM_VERSION "0.0.0"

#define IS_PAIR                       0x0001
#define ALL_SEGMENTS_MAPPED           0x0002
#define SEGMENT_UNMAPPED              0x0004
#define PAIRED_SEGMENT_UNMAPPED       0x0008
#define SEQ_REVERSED                  0x0010
#define NEXT_SEQ_REVERSED             0x0020
#define IS_FIRST_SEGMENT              0x0040
#define IS_LAST_SEGMENT               0x0080
#define IS_SECONDARY_ALIGNMENT        0x0100

#define UNDEF_NOT_PASSING_FILTERS     0x0200
#define UNDEF_OPTICAL_DUPLICATE       0x0400
#define UNDEF_SUPPLEMENTARY_ALIGNMENT 0x0800

#define UNDEF_FLAGS (UNDEF_NOT_PASSING_FILTERS | UNDEF_OPTICAL_DUPLICATE | UNDEF_SUPPLEMENTARY_ALIGNMENT)
#define UNUSED_FLAGS                  0xf000


#define SCORE_TO_MAPQ_FACTOR 1


void write_sam_alignment(genome_t* reference, int ref_seq_id, uint32_t position, char* qname, uint8_t mapq, char* cigar, int ref_id_next, int pnext, int tlen, char* seq, char* qual, uint16_t flag, int alignment_score, uint8_t template_independant_mapq);

void init_sam_file(char* sam_file_name, genome_t* reference);

void close_sam_file();

#endif // __SAMFILE_H__
