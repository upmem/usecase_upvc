/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#include "index.h"
#define _POSIX_C_SOURCE 200809L
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

#include "genome.h"
#include "parse_args.h"
#include "upvc.h"
#include "debug.h"
#include "profiling.h"

#define MAX_BUF_SIZE (1024)
#define GENOME_BINARY "genome.bin"
#define GENOME_VERSION (1)
#define GENOME_MAGIC (0x9e503e)

static genome_t genome = { .magic = GENOME_MAGIC, .version = GENOME_VERSION };

genome_t *genome_get() { return &genome; }

void genome_load()
{
    double start_time = my_clock();

    printf("%s:\n", __func__);
    char filename[FILENAME_MAX];
    char *index_folder = get_index_folder();
    sprintf(filename, "%s" GENOME_BINARY, index_folder);
    free(index_folder);
    FILE *f = fopen(filename, "r");
    CHECK_FILE(f, filename);

    xfer_file((uint8_t *)&genome, sizeof(genome_t), f, xfer_read);

    assert(genome.magic == GENOME_MAGIC
        && "Wrong header, make sure you have generated your MRAMs with the same version of UPVC that you are using.");
    assert(genome.version == GENOME_VERSION && "Could not load a genome generated with a different version of UPVC.");

    LOG_INFO("allocating genome data (%luMB + %luMB)\n", sizeof(int8_t) * genome.fasta_file_size/1000000, sizeof(int32_t) * genome.fasta_file_size/1000000);
    genome.data = (int8_t *)malloc(sizeof(int8_t) * genome.fasta_file_size);
    genome.mapping_coverage = (int32_t *)calloc(sizeof(int32_t), genome.fasta_file_size);
    assert(genome.data != NULL && genome.mapping_coverage != NULL);

    xfer_file((uint8_t *)genome.data, genome.fasta_file_size, f, xfer_read);

    fclose(f);

    printf("\t#seq: %d\n", genome.nb_seq);
    printf("\ttime: %lf s\n", my_clock() - start_time);
}

void genome_create()
{
    char *prefix = get_input_path();
    char genome_file_line[MAX_BUF_SIZE];

    double start_time = my_clock();

    printf("%s:\n", __func__);

    char filename[FILENAME_MAX];
    sprintf(filename, "%s.fasta", prefix);
    FILE *genome_file = fopen(filename, "r");
    CHECK_FILE(genome_file, filename);

    fseek(genome_file, 0, SEEK_END);
    genome.fasta_file_size = ftell(genome_file);
    rewind(genome_file);

    LOG_INFO("allocating genome data (%luMB)\n", sizeof(int8_t) * genome.fasta_file_size/1000000);
    genome.data = (int8_t *)malloc(sizeof(int8_t) * genome.fasta_file_size);
    assert(genome.data != NULL);
    genome.nb_seq = 0;
    genome.mapping_coverage = NULL;

    uint64_t current_data_idx = 0;
    while (fgets(genome_file_line, MAX_BUF_SIZE, genome_file) != NULL) {
        if (genome_file_line[0] == '>') { /* Commentary line with metadata */
            genome.pt_seq[genome.nb_seq] = current_data_idx;
            genome.len_seq[genome.nb_seq] = 0;
            genome_file_line[strlen(genome_file_line) - 1] = '\0';
            memcpy(genome.seq_name[genome.nb_seq], &genome_file_line[1],
                strlen(&genome_file_line[1]) > MAX_SEQ_NAME_SIZE ? MAX_SEQ_NAME_SIZE : strlen(&genome_file_line[1]));
            genome.nb_seq++;
        } else {
            for (unsigned int i = 0; i < strlen(genome_file_line) - 1; i++) {
                if (genome_file_line[i] != 'N') { /* A -> 0, C -> 1, G -> 3, T -> 2 */
                    genome.data[current_data_idx++] = (((int)genome_file_line[i]) >> 1) & 3;
                } else {
                    genome.data[current_data_idx++] = 4;
                }
                genome.len_seq[genome.nb_seq - 1]++;
            }
        }
    }
    fclose(genome_file);

    char *index_folder = get_index_folder();
    sprintf(filename, "%s" GENOME_BINARY, index_folder);
    free(index_folder);
    FILE *genome_binary = fopen(filename, "w");
    CHECK_FILE(genome_binary, filename);
    xfer_file((uint8_t *)&genome, sizeof(genome_t), genome_binary, xfer_write);
    xfer_file((uint8_t *)genome.data, genome.fasta_file_size, genome_binary, xfer_write);
    fclose(genome_binary);

    printf("\t#seq: %d\n", genome.nb_seq);
    printf("\ttime: %lf s\n", my_clock() - start_time);
}

void genome_free()
{
    free(genome.data);
    free(genome.mapping_coverage);
}

/**
 * Frequency table
 * 5 entries (A, C, T, G, -)
 * Each entry is a table of the size of the reference genome
 **/
static struct frequency_info* frequency_table[5];
static struct variants_codependence_info_list** variant_codependences;
#define NB_CODEPENDENCE_CHUNK_DIV (1<<10)
struct codependence_chunk* last_allocated_codependence_chunk[NB_CODEPENDENCE_CHUNK_DIV];
static bool init_frequency_table = false;
static bool init_codependence_table = false;

#define ALLOC_CHUNK_SIZE 100
struct codependence_chunk {
    struct codependence_chunk* previous_chunk;
    struct codependence_chunk* next_chunk;
    unsigned int next_slot_free;
    struct variants_codependence_info_list codependence_info_lists[ALLOC_CHUNK_SIZE];
};

static void allocate_new_codependence_chunk(int i) {
    STAT_RECORD_START(STAT_ALLOCATE_NEW_CHUNK);
    LOG_INFO("allocating codependence_chunk (%lu)\n", sizeof(struct codependence_chunk));
    struct codependence_chunk* new_chunk = calloc(1, sizeof(struct codependence_chunk));
    assert(new_chunk != NULL);
    new_chunk->previous_chunk = last_allocated_codependence_chunk[i];
    new_chunk->next_slot_free = 0;
    last_allocated_codependence_chunk[i] = new_chunk;
    STAT_RECORD_LAST_STEP(STAT_ALLOCATE_NEW_CHUNK, 0);
}

static struct variants_codependence_info_list* get_new_codependence_info_list(int i) {
    STAT_RECORD_START(STAT_GET_NEW_CODEPENDENCE_INFO);
    if (last_allocated_codependence_chunk[i]->next_slot_free >= ALLOC_CHUNK_SIZE) {
        allocate_new_codependence_chunk(i);
    }
    STAT_RECORD_LAST_STEP(STAT_GET_NEW_CODEPENDENCE_INFO, 0);
    return &(last_allocated_codependence_chunk[i]->codependence_info_lists[(last_allocated_codependence_chunk[i]->next_slot_free)++]);
}

void add_codependence_info(struct variants_codependence_info_list** next_variants_info_list, int16_t other_index_delta, uint8_t current_letter, uint8_t other_letter, unsigned int genome_size, pthread_mutex_t* mutex)
{
    STAT_RECORD_START(STAT_ADD_CODEPENDENCE_INFO);
    //struct variants_codependence_info_list** next_variants_info_list = &(freq_info->variants_codependence_list);
    struct variants_codependence_info_list* last_variants_info_list = NULL;
    int16_t key = (int16_t) ((other_index_delta<<4) | ((other_letter&0x3)<<2) | (current_letter&0x3));
    for (;*next_variants_info_list != NULL; next_variants_info_list = &((*next_variants_info_list)->next_list)) {
        for (int i=0; i<COD_LIST_SIZE; i++) {
            if ((*next_variants_info_list)->content[i].key == key) {
                if ((*next_variants_info_list)->content[i].codependence_count < UINT8_MAX)
                        (*next_variants_info_list)->content[i].codependence_count++;
                STAT_RECORD_LAST_STEP(STAT_ADD_CODEPENDENCE_INFO, 0);
                return;
            }
        }
        last_variants_info_list = *next_variants_info_list;
    }
    if (last_variants_info_list != NULL) {
        for (int i=0; i<COD_LIST_SIZE; i++) {
            if (last_variants_info_list->content[i].key == 0) {
                last_variants_info_list->content[i].key = key;
                last_variants_info_list->content[i].codependence_count = 1;
                STAT_RECORD_LAST_STEP(STAT_ADD_CODEPENDENCE_INFO, 1);
                return;
            }
        }
    }
    int chunk_div_index = (other_index_delta*NB_CODEPENDENCE_CHUNK_DIV)/genome_size;
    STAT_RECORD_STEP(STAT_ADD_CODEPENDENCE_INFO, 2);
    pthread_mutex_lock(mutex);
    STAT_RECORD_STEP(STAT_ADD_CODEPENDENCE_INFO, 3);
    struct variants_codependence_info_list* new_info_list = get_new_codependence_info_list(chunk_div_index);
    pthread_mutex_unlock(mutex);
    new_info_list->next_list = NULL;
    new_info_list->content[0].key = key;
    new_info_list->content[0].codependence_count = 1;
    *next_variants_info_list = new_info_list;
    STAT_RECORD_LAST_STEP(STAT_ADD_CODEPENDENCE_INFO, 4);
}

struct frequency_info** get_frequency_table() { 

  if(!init_frequency_table) {
    init_frequency_table = true;
    LOG_INFO("allocating frequency_table (5x%luMB)\n", sizeof(struct frequency_info)*genome.fasta_file_size/1000000);
    // allocate frequency_table on first call
    for(int i = 0; i < 5; ++i) {
      frequency_table[i] = (struct frequency_info*)calloc(genome.fasta_file_size, sizeof(struct frequency_info));
      assert(frequency_table[i] != NULL);
    }
  }
  return frequency_table; 
}

struct variants_codependence_info_list** get_codependence_table(pthread_mutex_t* mutex) {
    if(!init_codependence_table) {
        pthread_mutex_lock(mutex);
        if (!init_codependence_table) {
            LOG_INFO("allocating variant codependences lists (%luMB)\n", genome.fasta_file_size*sizeof(void*)/1000000);
            variant_codependences = (struct variants_codependence_info_list**)calloc(genome.fasta_file_size, sizeof(void*));
            assert(variant_codependences != NULL);
            for (int i=0; i<NB_CODEPENDENCE_CHUNK_DIV; i++) {
                last_allocated_codependence_chunk[i] = NULL;
                allocate_new_codependence_chunk(i);
            }
            init_codependence_table = true;
        }
        pthread_mutex_unlock(mutex);
    }
    return variant_codependences;
}

void free_frequency_table() {
  if(init_frequency_table) {
    for(int i = 0; i < 5; ++i) {
      free(frequency_table[i]);
    }
    init_frequency_table = false;
  }
}

void free_codependence_chunks() {
  // !!! This function does not remove pointers which are now pointing to unallocated regions
  LOG_INFO("deallocating codependence chunks\n");
  if(init_codependence_table) {
    for (int i=0; i<NB_CODEPENDENCE_CHUNK_DIV; i++) {
        struct codependence_chunk* to_delete_next = last_allocated_codependence_chunk[i];
        for (;last_allocated_codependence_chunk[i] != NULL; last_allocated_codependence_chunk[i] = to_delete_next) {
            to_delete_next = last_allocated_codependence_chunk[i]->previous_chunk;
            free(last_allocated_codependence_chunk[i]);
        }
    }
    init_codependence_table = false;
  }
}

