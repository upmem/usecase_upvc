#include "parse_args.h"
#include "samfile.h"

#include <assert.h>
#include <limits.h>
#include <pthread.h>
#include <stdio.h>

static FILE* sam_file = NULL;

static void write_sam_header(genome_t* reference) {
    fprintf(sam_file, "@HD\tVN:" SAM_VERSION "\n");

    for (uint32_t i=0; i<reference->nb_seq; i++) {
        // TODO: Maybe add AS tag (genome assembly identifier)
        fprintf(sam_file, "@SQ\tSN:%s\tLN:%lu\n", reference->seq_name[i], reference->len_seq[i]);
    }

    // TODO: Add CL tag (for command line used to call program)
    fprintf(sam_file, "@PG\tID:" PROGRAM_ID "\tPN:" PROGRAM_NAME "\tDS:" PROGRAM_DESCRIPTION "\tVN:" PROGRAM_VERSION "\tCL:%s\n", get_command_line());

    fprintf(sam_file, "@CO\tqname\tflag\trname\tpos\tmapq\tcigar\trnext\tpnext\ttlen\tseq\tqual\toptional_fields\n");
}

static char default_string[2] = "*";

static pthread_mutex_t sam_file_mutex = PTHREAD_MUTEX_INITIALIZER;

void write_sam_alignment(genome_t* reference, int ref_seq_id, uint32_t position, char* qname, uint8_t mapq, char* cigar, int ref_id_next, int pnext, int tlen, char* seq, char* qual, uint16_t flag, int alignment_score, uint8_t template_independant_mapq) {
    assert(sam_file != NULL);
    assert((flag & (UNUSED_FLAGS | UNDEF_FLAGS)) == 0);
    if (flag&IS_SECONDARY_ALIGNMENT) {
        // replace seq and qual with "*" for secondary alignments to reduce file size)
        seq = default_string;
        qual = default_string;
    }


    pthread_mutex_lock(&sam_file_mutex);
    fprintf(sam_file, "%s\t", qname);
    fprintf(sam_file, "%hu\t", flag);
    if (position == 0) {
        // If position unavailable, assume sequence unavailable too.
        fprintf(sam_file, "*\t");
    } else {
        fprintf(sam_file, "%s\t", reference->seq_name[ref_seq_id]);
    }
    fprintf(sam_file, "%u\t", position);
    fprintf(sam_file, "%hhu\t", mapq);
    fprintf(sam_file, "%s\t", cigar);
    if (pnext == 0) {
        // If pnext unavailable, assume seq_id unavailable too.
        fprintf(sam_file, "*\t");
    } else {
        if (ref_id_next == ref_seq_id) {
            fprintf(sam_file, "=\t");
        } else {
            fprintf(sam_file, "%s\t", reference->seq_name[ref_id_next]);
        }
    }
    fprintf(sam_file, "%u\t", pnext);
    fprintf(sam_file, "%i\t", tlen);
    fprintf(sam_file, "%s\t", seq);
    fprintf(sam_file, "%s\t", qual);
    
    // Optional fields
    if (alignment_score != INT_MAX) {
        fprintf(sam_file, "AS:i:%i\t", alignment_score);
    }
    if (template_independant_mapq != 255) {
        fprintf(sam_file, "SM:i:%i\t", template_independant_mapq);
    }
    fprintf(sam_file, "TC:i:2\n"); // This only works for read pairs for now.
    pthread_mutex_unlock(&sam_file_mutex);
}

void init_sam_file(char* sam_file_name, genome_t* reference) {
    sam_file = fopen(sam_file_name, "w");
    write_sam_header(reference);
}

void close_sam_file() {
    fclose(sam_file);
}
