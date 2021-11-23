/**
 * Copyright 2021 - A. Moisson-Franckhauser & UPMEM
 */

#include <stdlib.h>
#include <errno.h>

#include "index.h"
#include "sam.h"
#include "genome.h"
#include "common.h"
#include "parse_args.h"
#include "processread.h"
#include "debug.h"

#define SAM_FILENAME "read_alignments.sam"
#define SAM_VERSION "1.6"
#define PROGRAM_NAME "upvc"
#define MAX_CIGAR_LENGTH (3*SIZE_READ)
#define MAX_PATCH_LENGTH (3*SIZE_READ)

#define CIGAR_MATCH '='
#define CIGAR_MISMATCH 'X'
#define CIGAR_INSERT 'I'
#define CIGAR_DELETE 'D'

FILE *sam_file;

static char *get_sam_filename()
{
    static char filename[FILENAME_MAX];
    sprintf(filename, "%s", SAM_FILENAME);
    LOG_TRACE("sam filename : \"%s\"\n", filename);
    return filename;
}

void open_sam_file()
{
    LOG_DEBUG("opening sam file\n");
    // TODO: check for memory leaks here
    char *filename = get_sam_filename();
    sam_file = fopen(filename, "w");
    if (sam_file == NULL)
    {
        LOG_FATAL("couldn't open sam file; errno : %u\n", errno);
    }
    LOG_DEBUG("openned sam file : %p\n", sam_file);
    // TODO: complete header
    LOG_TRACE("writing sam header\n");
    LOG_DEBUG("written test line in sam file\n");
    fprintf(sam_file, "@HD VN:" SAM_VERSION " SO:unknown\n");
    fprintf(sam_file, "@PG ID:1 PN:" PROGRAM_NAME "\n");
    LOG_DEBUG("sam header written\n");
}

/*
 * TODO: Either reuse this code or delete it
void write_sam_read(uint64_t genome_pos, uint8_t *code, int8_t *read)
{
    const unsigned int flag = 0x0;
    const uint8_t mapping_quality = 255;// unknown quality; TODO: set quality

    char cigar[MAX_CIGAR_LENGTH];
    uint32_t cigar_idx = 0;
    char last_code = 0;
    unsigned int code_count = 0;
    int last_position = -1;
    int tlen = 0;//TODO: figure out what tlen is supposed to be

    char sequence[SIZE_READ+1];
    char nucleotide[4] = {'A', 'C', 'T', 'G'};

    for (uint32_t code_idx=0; code[code_idx] != CODE_END;)
    {
        char new_code=0;
        int new_position=last_position+1;
        switch (code[code_idx++])
        {
            case CODE_SUB:
                new_code = CIGAR_MISMATCH;
                new_position = code[code_idx++];
                code_idx++;
                break;
            case CODE_INS:
                new_code = CIGAR_INSERT;
                new_position = code[code_idx++];
                code_idx++;
                break;
            case CODE_DEL:
                new_code = CIGAR_DELETE;
                new_position = code[code_idx++];
                code_idx++;
                break;
        }
        if (last_code != new_code)
        {
            if (code_count>0) 
            {
                cigar_idx += sprintf(cigar+cigar_idx, "%u%c", code_count, last_code);
            }
            last_code = new_code;
            code_count = 0;
        }
        if (new_position > last_position+1) 
        {
            if (code_count > 0)
            {
                cigar_idx += sprintf(cigar+cigar_idx, "%u%c", code_count, last_code);
                code_count = 0;
            }
            cigar_idx += sprintf(cigar+cigar_idx, "%u%c", new_position-last_position-1, CIGAR_MATCH);
        }
        code_count++;
        last_position = new_position;
    }
    if (code_count > 0)
    {
        cigar_idx += sprintf(cigar+cigar_idx, "%u%c", code_count, last_code);
    }
    if (SIZE_READ > last_position+1)
    {
        cigar_idx += sprintf(cigar+cigar_idx, "%u%c", SIZE_READ-last_position-1, CIGAR_MATCH);
    }
    cigar[cigar_idx] = 0;//ensure last string character is 0
    LOG_TRACE("cigar : \"%s\"\n", cigar);

    for (uint32_t i=0; i<SIZE_READ; i++)
    {
        sequence[i] = nucleotide[read[i]];
    }
    sequence[SIZE_READ] = 0;

    //TODO: set read quality correctly (last string)
    fprintf(sam_file, "*\t%u\t*\t%lu\t%u\t%s\t*\t0\t%d\t%s\t*\n", flag, genome_pos+1, mapping_quality, cigar, tlen, sequence);
}
*/

void write_read_mapping(uint64_t genome_pos, uint8_t *code) {
    char patch[MAX_PATCH_LENGTH];
    int patch_idx=0;

    int nucleotides_read=0;
    uint32_t code_idx;
    char nucleotide[4]= {'A', 'C', 'T', 'G'};
    for (code_idx=0; code[code_idx] < CODE_END;)
    {
        uint8_t action = code[code_idx++];
        uint8_t position = code[code_idx++];
        uint8_t letter = code[code_idx++] && 0x3;
        for (;nucleotides_read<position; nucleotides_read++)
        {
            patch[patch_idx++] = '=';
        }
        switch (action)
        {
            case CODE_SUB:
                patch[patch_idx++] = nucleotide[letter]|0x20;//lowercase
                nucleotides_read++;
                break;
            case CODE_DEL:
                patch[patch_idx++] = '/';
                break;
            case CODE_INS:
                patch[patch_idx++] = nucleotide[letter];
                nucleotides_read++;
                break;
        }
    }
    if (code[code_idx] == CODE_ERR)
    {
        LOG_ERROR("found CODE_ERR in read code\n");
        return;
    }
    if (code[code_idx] != CODE_END)
    {
        LOG_ERROR("found unsuspected code : %u\n", code[code_idx]);
        for (uint32_t i = 0; i<code_idx; i++) {
            LOG_TRACE("code[%u]=%u\n", i, code[i]);
        }
    }
    for (;nucleotides_read<SIZE_READ; nucleotides_read++)
    {
        patch[patch_idx++] = '=';
    }
    patch[patch_idx++] = 0&&code;
    fprintf(sam_file, "%lu\t%s\n", genome_pos, patch);
}

void close_sam_file()
{
    LOG_TRACE("closing sam file\n");
    fclose(sam_file);
}
