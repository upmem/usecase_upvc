#include <stdio.h>
#include <string.h>

#include "parse_args.h"
#include "getread.h"
#include "genome.h"
#include "upvc.h"
#include "index.h"
#include "code.h"

int main(int argc, char *argv[])
{
        genome_t *ref_genome;
        times_ctx_t times_ctx;
        int *seed_counter_ref, *seed_counter_input;
        reads_info_t reads_info;
        int8_t *reads_buffer;
        int nb_read;
        char filename[1024];
        FILE *fipe1, *fipe2, *f_csv;
        char *input_prefix;
        unsigned long long nb_compute_total = 0ULL;
        unsigned long long nb_seed_ref = 0ULL;
        unsigned long long nb_seed_input = 0ULL;

        validate_args(argc, argv);

        input_prefix = get_input_path();
        reads_info.size_read = get_read_size(get_input_pe1());
        reads_info.size_neighbour_in_bytes = (reads_info.size_read - SIZE_SEED) / 4;
        ref_genome = get_genome(get_input_fasta(), &times_ctx);

        seed_counter_ref = (int *)calloc(NB_SEED, sizeof(int));
        seed_counter_input = (int *)calloc(NB_SEED, sizeof(int));
        for (int each_seq = 0; each_seq < ref_genome->nb_seq; each_seq++) {
                uint64_t sequence_start_idx = ref_genome->pt_seq[each_seq];
                for (uint64_t sequence_idx = 0;
                     sequence_idx < ref_genome->len_seq[each_seq] - reads_info.size_neighbour_in_bytes - SIZE_SEED + 1;
                     sequence_idx++) {
                        int seed_code = code_seed(&ref_genome->data[sequence_start_idx + sequence_idx]);
                        if (seed_code >= 0) {
                                seed_counter_ref[seed_code]++;
                        }
                }
        }

        sprintf(filename, "%s_PE1.fastq", input_prefix);
        fipe1 = fopen(filename, "r");
        sprintf(filename, "%s_PE2.fastq", input_prefix);
        fipe2 = fopen(filename, "r");

        reads_buffer = (int8_t *) malloc(sizeof(int8_t) * MAX_READS_BUFFER * reads_info.size_read);
        do {
                nb_read = get_reads(fipe1, fipe2, reads_buffer, &times_ctx, &reads_info);
                for (int each_read = 0; each_read < nb_read; each_read++) {
                        int8_t *curr_read = &reads_buffer[each_read * reads_info.size_read];
                        int seed_code = code_seed(curr_read);
                        if (seed_code >= 0) {
                                seed_counter_input[seed_code]++;
                                nb_compute_total += seed_counter_ref[seed_code];
                        }
                }
        } while (nb_read != 0);

        free(reads_buffer);
        fclose(fipe1);
        fclose(fipe2);

        printf("nb compute total = %e\n", (double)nb_compute_total);

        sprintf(filename, "%s_seed_counter.csv", input_prefix);
        f_csv = fopen(filename, "w");

        for (int each_seed = 0; each_seed < NB_SEED; each_seed++) {
                sprintf(filename, "%i, %i\n", seed_counter_ref[each_seed], seed_counter_input[each_seed]);
                fwrite(filename, sizeof(char), strlen(filename), f_csv);
                nb_seed_ref += seed_counter_ref[each_seed];
                nb_seed_input += seed_counter_input[each_seed];
        }

        printf("nb_seed_ref = %e\n"
               "nb_seed_input = %e\n",
               (double)nb_seed_ref,
               (double)nb_seed_input);

        fclose(f_csv);
        free(seed_counter_ref);
        free(seed_counter_input);
        free_args();

        return 0;
}
