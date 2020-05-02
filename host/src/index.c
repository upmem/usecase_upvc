/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#define _GNU_SOURCE
#include "index.h"
#include "genome.h"
#include "mram_dpu.h"
#include "parse_args.h"
#include "upvc.h"

#include <assert.h>
#include <errno.h>
#include <fcntl.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "common.h"

#define CODE_SIZE (4)
#define MASK (CODE_SIZE - 1)

int code_seed(int8_t *sequence)
{
    int seed = 0;
    for (int i = 0; i < SIZE_SEED; i++) {
        if (sequence[i] >= CODE_SIZE) {
            return -1;
        }
        seed = (seed * CODE_SIZE) + sequence[i];
    }
    return seed;
}

void code_neighbour(int8_t *sequence, int8_t *code)
{
    for (int i = 0; i < SIZE_NEIGHBOUR_IN_BYTES; i++) {
        int j = i * 4;
        code[i] = ((sequence[j + 3] & MASK) * (CODE_SIZE * CODE_SIZE * CODE_SIZE))
            + ((sequence[j + 2] & MASK) * (CODE_SIZE * CODE_SIZE)) + ((sequence[j + 1] & MASK) * (CODE_SIZE))
            + ((sequence[j] & MASK));
    }
}

void index_copy_neighbour(int8_t *dst, int8_t *src) { code_neighbour(&src[SIZE_SEED], dst); }

#define NB_SEED (1 << (SIZE_SEED << 1)) /* NB_SEED = 4 ^ (SIZE_SEED) */

#define MAX_SIZE_IDX_SEED (1000)

typedef struct hashtable_header {
    uint32_t magic;
    uint32_t version;
    uint32_t size_read;
    uint32_t size_seed;
    uint32_t nb_dpus;
    uint32_t unused;
    uint64_t nb_seed_total;
} hashtable_header_t;

#define INDEX_VERSION 1
static hashtable_header_t hashtable_header
    = { .magic = 0x1dec, .version = INDEX_VERSION, .size_read = SIZE_READ, .size_seed = SIZE_SEED };

static char *get_index_filename()
{
    static char filename[FILENAME_MAX];
    char *index_folder = get_index_folder();
    sprintf(filename, "%sindex.bin", index_folder);
    free(index_folder);
    return filename;
}

static unsigned int nb_indexed_dpu;
unsigned int index_get_nb_dpu() { return nb_indexed_dpu; }

static index_seed_t *index_seed;

index_seed_t *index_get(int8_t *read)
{
    index_seed_t *seed = &index_seed[code_seed(read)];
    if (seed->nb_nbr == 0 && seed->next == NULL)
        return NULL;
    else
        return seed;
}

typedef struct seed_counter {
    int nb_seed;
    int seed_code;
} seed_counter_t;

static int cmp_seed_counter(void const *a, void const *b)
{
    seed_counter_t *seed_counter_a = (seed_counter_t *)a;
    seed_counter_t *seed_counter_b = (seed_counter_t *)b;
    if (seed_counter_a->nb_seed < seed_counter_b->nb_seed) {
        return 1;
    } else {
        return -1;
    }
}

void index_load()
{
    double start_time = my_clock();
    printf("%s:\n", __func__);
    FILE *f = fopen(get_index_filename(), "r");
    CHECK_FILE(f, get_index_filename());

    hashtable_header_t header;
    fread(&header, sizeof(header), 1, f);
    assert(header.magic == hashtable_header.magic
        && "Wrong header, make sure you have generated your MRAMs with the same version of UPVC that you are "
           "using.");
    assert(header.version == hashtable_header.version && "Could not load an index generated with a different version of UPVC.");
    assert(header.size_read == hashtable_header.size_read && "Could not load an index generated with a different size of read.");
    assert(header.size_seed == hashtable_header.size_seed && "Could not load an index generated with a different size of seed.");

    index_seed = (index_seed_t *)malloc(header.nb_seed_total * sizeof(index_seed_t));
    assert(index_seed != NULL);

    xfer_file((uint8_t *)index_seed, sizeof(index_seed_t) * header.nb_seed_total, f, xfer_read);

    for (unsigned int each_seed = 0; each_seed < header.nb_seed_total; each_seed++) {
        if ((uintptr_t)index_seed[each_seed].next == UINTPTR_MAX) {
            index_seed[each_seed].next = NULL;
        } else {
            index_seed[each_seed].next = (struct index_seed *)((uintptr_t)index_seed[each_seed].next + (uintptr_t)index_seed);
        }
    }

    nb_indexed_dpu = header.nb_dpus;
    printf("\tnb_dpu: %u\n"
           "\tsize_read: %u\n"
           "\tsize_seed: %u\n",
        nb_indexed_dpu, header.size_read, header.size_seed);

    fclose(f);
    printf("\ttime: %lf s\n", my_clock() - start_time);
}

static int compute_nb_index_needed(int nb_seed) { return (nb_seed + MAX_SIZE_IDX_SEED - 1) / MAX_SIZE_IDX_SEED; }

#define INDEX_THREAD (16)
#define INDEX_THREAD_SLAVE (INDEX_THREAD - 1)
static seed_counter_t *seed_counter;
static pthread_barrier_t barrier;
static uint64_t nb_seed_total;

static void set_seed_counter(int thread_id)
{
    pthread_barrier_wait(&barrier);
    for (int i = thread_id; i < NB_SEED; i += INDEX_THREAD) {
        seed_counter[i].nb_seed = 0;
        seed_counter[i].seed_code = i;
    }

    static genome_t *ref_genome;
    if (thread_id == 0)
        ref_genome = genome_get();
    pthread_barrier_wait(&barrier);

    for (uint32_t i = 0; i < ref_genome->nb_seq; i++) {
        uint64_t sequence_start_idx = ref_genome->pt_seq[i];
        for (uint64_t sequence_idx = thread_id; sequence_idx < ref_genome->len_seq[i] - SIZE_NEIGHBOUR_IN_BYTES - SIZE_SEED + 1;
             sequence_idx += INDEX_THREAD) {
            int seed_code = code_seed(&ref_genome->data[sequence_start_idx + sequence_idx]);
            if (seed_code >= 0) {
                __sync_fetch_and_add(&seed_counter[seed_code].nb_seed, 1);
            }
        }
    }
    pthread_barrier_wait(&barrier);
}

static void init_index_seed(int thread_id)
{
    static uint64_t seed_offset = NB_SEED;

    pthread_barrier_wait(&barrier);
    for (int i = thread_id; i < NB_SEED; i += INDEX_THREAD) {
        if (seed_counter[i].nb_seed == 0) {
            index_seed[i].nb_nbr = 0;
            index_seed[i].next = NULL;
            continue;
        }

        int nb_index_needed = compute_nb_index_needed(seed_counter[i].nb_seed);
        int nb_neighbour_per_index = (seed_counter[i].nb_seed + nb_index_needed - 1) / nb_index_needed;
        index_seed[i].nb_nbr = seed_counter[i].nb_seed - ((nb_index_needed - 1) * nb_neighbour_per_index);
        index_seed[i].next = NULL;
        for (int j = 1; j < nb_index_needed; j++) {
            index_seed_t *seed = &index_seed[__sync_fetch_and_add(&seed_offset, 1)];
            seed->nb_nbr = nb_neighbour_per_index;
            seed->next = index_seed[i].next;
            index_seed[i].next = seed;
        }
    }
    pthread_barrier_wait(&barrier);
}

static void adjust_index_seed_next(int thread_id)
{
    pthread_barrier_wait(&barrier);
    for (unsigned int i = thread_id; i < nb_seed_total; i += INDEX_THREAD) {
        if (index_seed[i].next == NULL) {
            index_seed[i].next = (index_seed_t *)(UINTPTR_MAX);
        } else {
            index_seed[i].next = (index_seed_t *)((uintptr_t)index_seed[i].next - (uintptr_t)index_seed);
        }
    }
    pthread_barrier_wait(&barrier);
}
static void write_data(int thread_id)
{
    static genome_t *ref_genome;
    static volatile uint32_t seq_number_shared[INDEX_THREAD_SLAVE] = { 0 };
    static volatile uint64_t sequence_idx_shared[INDEX_THREAD_SLAVE] = { 0 };
    coords_and_nbr_t buffer;
    if (thread_id == 0)
        ref_genome = genome_get();
    pthread_barrier_wait(&barrier);
    if (thread_id == INDEX_THREAD_SLAVE) {
        double start_time = my_clock();
        while (true) {
            uint32_t slowest_thread_seq = UINT_MAX;
            uint64_t slowest_thread_sequence = ULONG_MAX;
            for (unsigned int each_thread = 0; each_thread < INDEX_THREAD_SLAVE; each_thread++) {
                uint32_t curr_seq = seq_number_shared[each_thread];
                uint32_t sequence_curr = sequence_idx_shared[each_thread];
                if (curr_seq < slowest_thread_seq) {
                    slowest_thread_seq = curr_seq;
                    slowest_thread_sequence = sequence_curr;
                } else if (curr_seq == slowest_thread_seq && sequence_curr < slowest_thread_sequence) {
                    slowest_thread_sequence = sequence_curr;
                }
            }
            double time = my_clock() - start_time;
            char time_str[FILENAME_MAX];
            if (time < 60.0) {
                sprintf(time_str, " - %us ", (unsigned int)time);
            } else {
                sprintf(time_str, " - %umin%us ", (unsigned int)time / 60, (unsigned int)time % 60);
            }
            if (slowest_thread_seq == ref_genome->nb_seq) {
                printf("\r\t\t100.00%%#%u%s\n", ref_genome->nb_seq, time_str);
                break;
            }
            printf("\r\t\t%2.2f%%#%u%s", ((float)slowest_thread_sequence) * 100.0 / ref_genome->len_seq[slowest_thread_seq],
                slowest_thread_seq, time_str);
            fflush(stdout);
            usleep(1000000);
        }
    } else {
        for (uint32_t seq_number = 0; seq_number < ref_genome->nb_seq; seq_number++, seq_number_shared[thread_id] = seq_number) {
            uint64_t sequence_start_idx = ref_genome->pt_seq[seq_number];
            for (uint64_t sequence_idx = thread_id;
                 sequence_idx < ref_genome->len_seq[seq_number] - SIZE_NEIGHBOUR_IN_BYTES - SIZE_SEED + 1;
                 sequence_idx += INDEX_THREAD_SLAVE, sequence_idx_shared[thread_id] = sequence_idx) {
                index_seed_t *seed;
                int align_idx;
                int seed_code = code_seed(&ref_genome->data[sequence_start_idx + sequence_idx]);

                if (seed_code < 0) {
                    continue;
                }

                seed = &index_seed[seed_code];

                int total_nb_neighbour = 0;
                int32_t nb_seed = __sync_fetch_and_add(&seed_counter[seed_code].nb_seed, 1);
                while (seed != NULL) {
                    if (nb_seed < (int)seed->nb_nbr + total_nb_neighbour)
                        break;
                    total_nb_neighbour += seed->nb_nbr;
                    seed = seed->next;
                }
                align_idx = seed->offset + nb_seed - total_nb_neighbour;

                buffer.coord.seq_nr = seq_number;
                buffer.coord.seed_nr = sequence_idx;
                code_neighbour(&ref_genome->data[sequence_start_idx + sequence_idx + SIZE_SEED], (int8_t *)&buffer.nbr);
                write_vmi(seed->num_dpu, align_idx, &buffer);
            }
        }
    }
    pthread_barrier_wait(&barrier);
}

static void *index_create_slave_fct(void *args)
{
    uint32_t thread_id = (uint32_t)(uintptr_t)args;

    set_seed_counter(thread_id);
    init_index_seed(thread_id);
    write_data(thread_id);
    adjust_index_seed_next(thread_id);

    return NULL;
}

void index_create_folder()
{
    char *index_folder = get_index_folder();
    if (access(index_folder, F_OK) == 0) {
        ERROR_EXIT(ERR_INDEX_FOLDER_ALREAD_EXIST, "'%s' already exits, please save it elsewhere or erase it", index_folder);
    }
    if (mkdir(index_folder, 0755) != 0) {
        ERROR_EXIT(ERR_INDEX_FOLDER_MKDIR_FAILED,
            "Could not create '%s', make sure you have the permission to create folder and file here", index_folder);
    }
    free(index_folder);
}

void index_create()
{
    unsigned int nb_dpu = get_nb_dpu();
    double start_time = my_clock();
    printf("%s(%i):\n", __func__, nb_dpu);
    pthread_t thread_id[INDEX_THREAD_SLAVE];
    distribute_index_t *distribute_index_table;

    assert(pthread_barrier_init(&barrier, NULL, INDEX_THREAD) == 0);
    for (unsigned int each_thread = 0; each_thread < INDEX_THREAD_SLAVE; each_thread++) {
        assert(pthread_create(&thread_id[each_thread], NULL, index_create_slave_fct, (void *)(uintptr_t)each_thread) == 0);
    }

    {
        double init_seed_counter_time = my_clock();
        printf("\tInitialize the seed_counter table\n");
        seed_counter = (seed_counter_t *)malloc(NB_SEED * sizeof(seed_counter_t));
        set_seed_counter(INDEX_THREAD_SLAVE);

        printf("\t\ttime: %lf s\n", my_clock() - init_seed_counter_time);
    }

    {
        double alloc_index_seed_time = my_clock();
        printf("\tAllocating the index table\n");
        nb_seed_total = NB_SEED;
        for (int i = 0; i < NB_SEED; i++) {
            nb_seed_total += compute_nb_index_needed(seed_counter[i].nb_seed);
        }
        index_seed = (index_seed_t *)malloc(sizeof(index_seed_t) * nb_seed_total);
        assert(index_seed != NULL);
        printf("\t\tnb_seed_total=%lu\n"
               "\t\ttime: %lf s\n",
            nb_seed_total, my_clock() - alloc_index_seed_time);
    }

    {
        double create_init_link_all_seed_time = my_clock();
        printf("\tCreate, initialize and link together all the seed\n");
        init_index_seed(INDEX_THREAD_SLAVE);
        printf("\t\ttime: %lf s\n", my_clock() - create_init_link_all_seed_time);
    }

    {
        double sort_time = my_clock();
        printf("\tSort seed counter\n");
        qsort(seed_counter, NB_SEED, sizeof(seed_counter_t), cmp_seed_counter);
        printf("\t\ttime: %lf s\n", my_clock() - sort_time);
    }

    {
        double distribute_index_time = my_clock();
        printf("\tDistribute indexs between DPUs\n");

        distribute_index_table = (distribute_index_t *)calloc(nb_dpu, sizeof(distribute_index_t));
        assert(distribute_index_table != NULL);
        struct distribute_index_list head = TAILQ_HEAD_INITIALIZER(head);

        for (unsigned int i = 0; i < nb_dpu; i++) {
            distribute_index_table[i].dpu_id = i;
            TAILQ_INSERT_TAIL(&head, &distribute_index_table[i], entries);
        }

        for (int i = 0; i < NB_SEED; i++) {
            int seed_code = seed_counter[i].seed_code;
            int nb_seed_counted = seed_counter[i].nb_seed;
            index_seed_t *seed = &index_seed[seed_code];
            if (seed->nb_nbr == 0 && seed->next == NULL) {
                continue;
            }

            while (seed != NULL) {
                distribute_index_t *dpu = TAILQ_LAST(&head, distribute_index_list);
                TAILQ_REMOVE(&head, dpu, entries);
                seed->offset = dpu->size;
                seed->num_dpu = dpu->dpu_id;
                dpu->size += seed->nb_nbr;
                dpu->workload += (uint64_t)seed->nb_nbr * (uint64_t)nb_seed_counted;
                seed = seed->next;

                distribute_index_t *dpu_cmp;
                TAILQ_FOREACH(dpu_cmp, &head, entries)
                {
                    if (dpu->workload >= dpu_cmp->workload) {
                        TAILQ_INSERT_BEFORE(dpu_cmp, dpu, entries);
                        break;
                    }
                }
            }
        }

        printf("\t\ttime: %lf s\n", my_clock() - distribute_index_time);
    }

    /* Writing data in DPUs memories */
    {
        double write_in_memories_time = my_clock();
        printf("\tWriting data in DPUs memories\n");

        memset(seed_counter, 0, sizeof(seed_counter_t) * NB_SEED);

        init_vmis(nb_dpu, distribute_index_table);
        write_data(INDEX_THREAD_SLAVE);
        free_vmis(nb_dpu);

        free(seed_counter);
        free(distribute_index_table);
        printf("\t\ttime: %lf s\n", my_clock() - write_in_memories_time);
    }

    {
        double start_time = my_clock();
        printf("\tSaving index on disk\n");
        FILE *f = fopen(get_index_filename(), "w");
        CHECK_FILE(f, get_index_filename());

        hashtable_header.nb_seed_total = nb_seed_total;
        hashtable_header.nb_dpus = nb_dpu;
        fwrite(&hashtable_header, sizeof(hashtable_header_t), 1, f);

        adjust_index_seed_next(INDEX_THREAD_SLAVE);

        xfer_file((uint8_t *)index_seed, sizeof(index_seed_t) * nb_seed_total, f, xfer_write);

        fclose(f);
        printf("\t\ttime: %lf s\n", my_clock() - start_time);
    }

    for (unsigned int each_thread = 0; each_thread < INDEX_THREAD_SLAVE; each_thread++) {
        assert(pthread_join(thread_id[each_thread], NULL) == 0);
    }
    assert(pthread_barrier_destroy(&barrier) == 0);

    printf("\ttime: %lf s\n", my_clock() - start_time);
}

void index_free() { free(index_seed); }

char *get_index_folder()
{
    char *filename;
    assert(asprintf(&filename, "%s_index/", get_input_path()) > 0);
    assert(filename != NULL);
    return filename;
}
