/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#ifndef __GENOME_H__
#define __GENOME_H__

#include <stdint.h>
#include <pthread.h>

#define MAX_SEQ_GEN (24) // max number of chromosomes
#define MAX_SEQ_NAME_SIZE (8)

typedef struct {
    uint32_t magic;
    uint32_t version;
    uint32_t nb_seq; // nb of chromosomes (24)
    uint64_t pt_seq[MAX_SEQ_GEN]; // offset of each chromosome in data
    uint64_t len_seq[MAX_SEQ_GEN]; // length of chromosome in data
    uint64_t fasta_file_size;
    char seq_name[MAX_SEQ_GEN][MAX_SEQ_NAME_SIZE];
    int8_t *data; //genome de reference 1B = 1 nucleotide
    int32_t *mapping_coverage;
} genome_t;

void genome_create();

void genome_load();

void genome_free();

genome_t *genome_get();

struct frequency_info {

  float freq;
  unsigned int score;
  unsigned int unsure_score;
};

#pragma pack(push,1)
struct variants_codependence_info {
    int16_t key;
    uint8_t codependence_count;
};

#define COD_LIST_SIZE 4
struct variants_codependence_info_list {
    struct variants_codependence_info_list* next_list;
    struct variants_codependence_info content[COD_LIST_SIZE];
};
#pragma pack(pop)

void add_codependence_info(struct variants_codependence_info_list** next_variants_info_list, int16_t other_index_delta, uint8_t current_letter, uint8_t other_letter, unsigned int genome_size, pthread_mutex_t* mutex);

struct frequency_info** get_frequency_table(); 
struct variants_codependence_info_list** get_codependence_table(pthread_mutex_t* mutex);
void free_frequency_table();
void free_codependence_chunks();

#endif /* __GENOME_H__ */
