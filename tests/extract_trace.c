#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NB_THREAD (16)
#define MAX_STRING_SIZE (512)
#define MAX_BLOB_SIZE (1000000)
#define PC_FORMAT ("  PC=0x%x(%u) TID=%u CYCLE=%u\n")
#define EXPECTED_MATCH_PC (2)
#define CALL_FORMAT ("\tcall 0x%x\n")
#define EXPECTED_MATCH_CALL (1)
#define RET_FORMAT ("\tret 0x%x\n")
#define EXPECTED_MATCH_RET (1)

static bool goto_thread_info(char *str, unsigned int thread_id, FILE *input_fp)
{
    while (true) {
        int nb_match;
        unsigned int pc_addr;
        unsigned int tid;
        __attribute__((unused)) unsigned int unused1, unused2;
        if (fgets(str, MAX_STRING_SIZE, input_fp) == NULL) {
            return false;
        }
        nb_match = sscanf(str, PC_FORMAT, &pc_addr, &unused1, &tid, &unused2);
        if (nb_match == EXPECTED_MATCH_PC && tid == thread_id) {
            return true;
        }
    }
}

static bool save_thread_info(char *str, unsigned int call_addr_to_track, unsigned int thread_id, unsigned int *contain_call,
    FILE *input_fp, FILE *output_fp)
{
    char blob[MAX_BLOB_SIZE];
    unsigned int blob_size = 0;

    memcpy(&blob[blob_size], str, strlen(str));
    blob_size += strlen(str);
    assert(blob_size < MAX_BLOB_SIZE);

    while (true) {
        int nb_match;
        unsigned int pc_addr;
        unsigned int tid;
        __attribute__((unused)) unsigned int unused1, unused2;
        if (fgets(str, MAX_STRING_SIZE, input_fp) == NULL) {
            if (*contain_call) {
                fwrite(blob, sizeof(char), blob_size, output_fp);
            }
            return false;
        }

        nb_match = sscanf(str, PC_FORMAT, &pc_addr, &unused1, &tid, &unused2);

        if (nb_match == EXPECTED_MATCH_PC && tid != thread_id) {
            if (*contain_call) {
                fwrite(blob, sizeof(char), blob_size, output_fp);
            }
            return true;
        } else if (nb_match == EXPECTED_MATCH_PC && tid == thread_id) {
            if (!*contain_call) {
                blob_size = 0;
            }
            memcpy(&blob[blob_size], str, strlen(str));
            blob_size += strlen(str);
            assert(blob_size < MAX_BLOB_SIZE);
        } else if (nb_match != EXPECTED_MATCH_PC) {
            unsigned int addr;

            memcpy(&blob[blob_size], str, strlen(str));
            blob_size += strlen(str);
            assert(blob_size < MAX_BLOB_SIZE);

            nb_match = sscanf(str, CALL_FORMAT, &addr);
            if (nb_match == EXPECTED_MATCH_CALL && addr == call_addr_to_track) {
                if (*contain_call == 0) {
                    fprintf(output_fp, "################################\n");
                }
                *contain_call = *contain_call + 1;
                continue;
            } else if (nb_match == EXPECTED_MATCH_CALL && *contain_call > 0) {
                *contain_call = *contain_call + 1;
                continue;
            }

            nb_match = sscanf(str, RET_FORMAT, &addr);
            if (nb_match == EXPECTED_MATCH_RET) {
                if (*contain_call <= 0) {
                    continue;
                }
                *contain_call = *contain_call - 1;
                if (!*contain_call) {
                    fwrite(blob, sizeof(char), blob_size, output_fp);
                    fprintf(output_fp, "################################\n");
                    blob_size = 0;
                }
                continue;
            }
        }
    }
}

static void extract_thread(unsigned int call_addr_to_track, unsigned int thread_id, FILE *input_fp, FILE *output_fp)
{
    char str[MAX_STRING_SIZE];
    unsigned contain_call = 0;
    while (true) {
        if (!goto_thread_info(str, thread_id, input_fp)) {
            break;
        }
        if (!save_thread_info(str, call_addr_to_track, thread_id, &contain_call, input_fp, output_fp)) {
            break;
        }
    }
}

int main(__attribute__((unused)) int argc, char **argv)
{
    for (unsigned int i = 0; i < NB_THREAD; i++) {
        char output_file[MAX_STRING_SIZE];
        FILE *input_fp = fopen(argv[1], "r");

        printf("Extracting thread %u\n", i);

        sprintf(output_file, "%s.%u", argv[1], i);
        FILE *output_fp = fopen(output_file, "w");

        extract_thread(atoi(argv[2]), i, input_fp, output_fp);

        fclose(output_fp);
        fclose(input_fp);
    }
    return 0;
}
