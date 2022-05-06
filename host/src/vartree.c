/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#include <math.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "common.h"
#include "genome.h"
#include "parse_args.h"
#include "upvc.h"
#include "vartree.h"
#include "debug.h"


#define DUMP_FREQUENCY_TABLE false
#define FREQ_TABLE_DUMP_FILE_NAME "frequency_table.bin"
#if DUMP_FREQUENCY_TABLE
FILE *freq_table_dump_file;

void open_freq_table()
{
    static char filename[FILENAME_MAX];
    sprintf(filename, "%s", FREQ_TABLE_DUMP_FILE_NAME);
    freq_table_dump_file = fopen(filename, "w");
    if (freq_table_dump_file == NULL) {
        LOG_FATAL("couldn't open frequency table dumping file; errno: %u\n", errno);
    }
}

void close_freq_table()
{
    fclose(freq_table_dump_file);
}

//The pragma pack shouldn't have any effect here as of now but it might have an effect if struct content is modified
//In which case it should probably be kept.
#pragma pack(push,1)
struct freq_table_dump_entry_t{
    float freqs[5];
    uint16_t  scores[5];
};
#pragma pack(pop)

// freq_table_dump format:
//  0:            uint8_t                                 number of sequences = {n}
//  1 - n*4:      uint16_t[n]                             sequence start addresses
//  n*4+1 - end:  struct freq_table_dump_entry_t[...]     table entries
void write_freq_table_header(genome_t *ref_genome)
{
    uint8_t number_sequences = ref_genome->nb_seq;
    uint64_t sequence_offsets[256];
    uint64_t header_size = 1 + sizeof(uint64_t)* (uint64_t) number_sequences;
    uint64_t total_seq_length=0;
    for(uint8_t seq_number=0; seq_number<ref_genome->nb_seq; seq_number++) {
        sequence_offsets[seq_number] = header_size + total_seq_length*sizeof(struct freq_table_dump_entry_t);
        total_seq_length += ref_genome->len_seq[seq_number];
    }
    fwrite(&number_sequences, 1, 1, freq_table_dump_file);
    fwrite(sequence_offsets, sizeof(uint64_t), number_sequences, freq_table_dump_file);
}

void dump_freq_table_entry(int address, struct frequency_info **frequency_table)
{
    struct freq_table_dump_entry_t entry;
    for (int i=0; i<5; i++) {
        entry.freqs[i] = frequency_table[i][address].freq;
        entry.scores[i] = frequency_table[i][address].score;
    }
    fwrite(&entry, sizeof(struct freq_table_dump_entry_t), 1, freq_table_dump_file);
}

void dump_freq_table(genome_t *ref_genome, struct frequency_info **frequency_table)
{
    open_freq_table();
    write_freq_table_header(ref_genome);
    for (uint32_t seq_number = 0; seq_number < ref_genome->nb_seq; seq_number++) {
        for (uint64_t seq_position = 0; seq_position < ref_genome->len_seq[seq_number]; seq_position++) {
            dump_freq_table_entry(ref_genome->pt_seq[seq_number] + seq_position, frequency_table);
        }
    }
    close_freq_table();
}

#endif // DUMP_FREQUENCY_TABLE

static variant_t **variant_list[MAX_SEQ_GEN] = { NULL };
static pthread_mutex_t mutex;

void variant_tree_init()
{
    genome_t *genome = genome_get();
    pthread_mutex_init(&mutex, NULL);
    for (unsigned int each_seq = 0; each_seq < genome->nb_seq; each_seq++) {
        LOG_INFO("allocating variant_list (%luMB)\n", sizeof(variant_t*) * genome->len_seq[each_seq]);
        variant_list[each_seq] = (variant_t **)calloc(genome->len_seq[each_seq], sizeof(variant_t *));
    }
}

void variant_tree_insert(variant_t *var, uint32_t seq_nr, uint32_t offset_in_chr)
{
    pthread_mutex_lock(&mutex);
    variant_t **entry = &variant_list[seq_nr][offset_in_chr];
    variant_t *vars = *entry;
    while (vars != NULL) {
        if (!strcmp(vars->ref, var->ref) && !strcmp(vars->alt, var->alt)) {
            vars->depth++;
            vars->score += var->score;
            free(var);
            goto end;
        }
        vars = vars->next;
    }
    var->next = *entry;
    *entry = var;

end:
    pthread_mutex_unlock(&mutex);
}

void variant_tree_free()
{
    genome_t *genome = genome_get();
    pthread_mutex_destroy(&mutex);
    for (unsigned int each_seq = 0; each_seq < genome->nb_seq; each_seq++) {
        for (unsigned int i = 0; i < genome->len_seq[each_seq]; i++) {
            variant_t *tmp = variant_list[each_seq][i];
            while (tmp != NULL) {
                variant_t *to_free = tmp;
                tmp = tmp->next;
                free(to_free);
            }
        }
        free(variant_list[each_seq]);
    }
}

typedef struct {
    uint32_t percentage;
    uint32_t score;
} depth_filter_t;

#if 0
#if (SIZE_READ == 120)
depth_filter_t sub_filter[] = {
    [3] = { 15, 16 },
    [4] = { 17, 17 },
    [5] = { 18, 18 },
    [6] = { 20, 18 },
    [7] = { 21, 20 },
    [8] = { 22, 21 },
    [9] = { 22, 21 },
    [10] = { 24, 21 },
    [11] = { 24, 21 },
    [12] = { 28, 21 },
    [13] = { 29, 22 },
    [14] = { 29, 23 },
    [15] = { 32, 24 },
    [16] = { 32, 25 },
    [17] = { 35, 25 },
    [18] = { 35, 25 },
    [19] = { 35, 25 },
    [20] = { 40, 25 },
};

depth_filter_t indel_filter[] = {
    [2] = { 10, 16 },
    [3] = { 12, 21 },
    [4] = { 13, 21 },
    [5] = { 14, 22 },
    [6] = { 14, 22 },
    [7] = { 1, 23 },
    [8] = { 1, 25 },
    [9] = { 1, 25 },
    [10] = { 1, 30 },
    [11] = { 1, 40 },
};
#elif (SIZE_READ == 150) || (SIZE_READ==148)
depth_filter_t sub_filter[] = {
    [3] = { 15, 16 },
    [4] = { 17, 20 },
    [5] = { 18, 20 },
    [6] = { 20, 21 },
    [7] = { 21, 21 },
    [8] = { 22, 21 },
    [9] = { 24, 22 },
    [10] = { 25, 23 },
    [11] = { 27, 23 },
    [12] = { 27, 25 },
    [13] = { 29, 25 },
    [14] = { 30, 27 },
    [15] = { 31, 27 },
    [16] = { 34, 27 },
    [17] = { 34, 27 },
    [18] = { 34, 29 },
    [19] = { 35, 29 },
    [20] = { 40, 29 },
};

depth_filter_t indel_filter[] = {
    [2] = { 9, 21 },
    [3] = { 12, 22 },
    [4] = { 12, 22 },
    [5] = { 13, 24 },
    [6] = { 15, 25 },
    [7] = { 17, 25 },
    [8] = { 18, 25 },
    [9] = { 2, 26 },
    [10] = { 1, 27 },
    [11] = { 1, 40 },
};
#else
#error "Filter not defined for this size of read"
#endif
#endif

#if 0
static bool homopolymer(int8_t *seq, int offset)
{
    for (int i = 0; i < offset - 1; i++) {
        if (seq[i] != seq[i + 1]) {
            return false;
        }
    }
    return true;
}
#endif

static bool print_variant_tree(variant_t *var, uint32_t seq_nr, uint64_t seq_pos, genome_t *ref_genome, FILE *vcf_file)
{
    char *chr = ref_genome->seq_name[seq_nr];
    uint64_t genome_pos = ref_genome->pt_seq[seq_nr] + seq_pos;
    uint32_t cov = ref_genome->mapping_coverage[genome_pos];
    uint32_t depth = var->depth;
    uint32_t score = var->score / depth;
    // Note: commenting out the old version of variant calling, now using the frequency table
    //uint32_t percentage = 100;
    //if (cov != 0) {
    //    percentage = depth * 100 / cov;
    //}

    uint32_t ref_len = strlen(var->ref);
    uint32_t alt_len = strlen(var->alt);
    //if (ref_len > alt_len && percentage <= 25 && homopolymer(&ref_genome->data[genome_pos - 12], 12)) {
    //    return false;
    //}

    if (get_no_filter())
        goto print;

    if (ref_len == alt_len) { /* SUBSTITUTION */
        //if (depth < 3) {
        //    return false;
        //} else if (depth > 20) {
        //    depth = 20;
        //}
        //if (!(score <= sub_filter[depth].score && percentage >= sub_filter[depth].percentage)) {
        //    return false;
        //}
        if (depth > 20) {
            depth = 20;
        }
    } else { /* INSERTION OR DELETION */
        if (depth < 2) {
            return false;
        } else if (depth > 11) {
            depth = 11;
        }
        //if (!(score <= indel_filter[depth].score && percentage >= indel_filter[depth].percentage)) {
        //    return false;
        //}
    }

print:
    //TODO
    fprintf(vcf_file, "%s\t%lu\t.\t%s\t%s\t.\t.\tDEPTH=%d;COV=%d;SCORE=%d\n", chr, seq_pos+1, var->ref, var->alt, var->depth, cov,
        score);

    return true;
}

/**
  Few configurations suggested by Bertil
  (i) D=3: 25%, D=4: 20%, D=5: 15%; D>=6: 10%
  (ii) D=2: 30%, D=3: 25%, D=4: 20%, D=5: 15%; D>=6: 10%
  (iii) D=3: 20%, D=4: 15%, D>=5: 10%;
 **/

__attribute__((unused)) uint32_t depth_filter1(float freq) {
  if(freq < 10.0f)
    return UINT_MAX;
  if(freq < 15.0f)
    return 6;
  if(freq < 20.0f)
    return 5;
  if(freq < 25.0f)
    return 4;
  return 3;
}

__attribute__((unused)) uint32_t depth_filter2(float freq) {
  if(freq < 10.0f)
    return UINT_MAX;
  if(freq < 15.0f)
    return 6;
  if(freq < 20.0f)
    return 5;
  if(freq < 25.0f)
    return 4;
  if(freq < 30.0f)
    return 3;
  return 2;
}

__attribute__((unused)) uint32_t depth_filter3(float freq) {
  if(freq < 10.0f)
    return UINT_MAX;
  if(freq < 15.0f)
    return 5;
  if(freq < 20.0f)
    return 4;
  return 3;
}

__attribute__((unused)) uint32_t depth_filter_a(float freq) {
    if(freq >= 20)
        return 3;
    if(freq >= 15)
        return 5;
    return UINT_MAX;
}

__attribute__((unused)) uint32_t depth_filter_permissive(float freq) {
    if (freq< 10.0f) {
        return 3;
    }
    if (freq<20.0f) {
        return 2;
    }
    return 1;
}

__attribute__((unused)) uint32_t depth_filter_fixed_3(float freq) {

  if(freq < 20.0f)
    return UINT_MAX;
  return 3;
}

__attribute__((unused)) uint32_t depth_filter_fixed_3_f15(float freq) {

  if(freq < 15.0f)
    return UINT_MAX;
  return 3;
}

#define AFFINE_B 1.98142
#define AFFINE_A 0.164761

float reverse_filter(uint32_t score) {
    return AFFINE_A*(float)score + AFFINE_B;
}

FILE * dbg_file = NULL;
FILE * sub_file = NULL;

static void get_most_frequent_variant(genome_t * ref_genome, struct frequency_info ** frequency_table, uint32_t seq_number, uint64_t seq_position, variant_t * results) {

  static char nucleotide[4] = { 'A', 'C', 'T', 'G' };// FIXME : const

  uint64_t genome_pos = ref_genome->pt_seq[seq_number] + seq_position;

  float total = 0;
  for(int i = 0; i < 5; ++i) {
    total += frequency_table[i][genome_pos].freq; 
    results[i].depth = 0;
    results[i].score = 0;
  }
  if(total == 0) 
    return ;//results;

  for(int i = 0; i < 5; ++i) {
    float freq = frequency_table[i][genome_pos].freq;
    //uint32_t score = frequency_table[i][genome_pos].score;
    if(i == ref_genome->data[genome_pos]) {
        continue; // not a variant if the same nucleotide as in reference genome
    }

    // if frequency and depth pass the threshold, consider it a variant
    if(freq > reverse_filter(total)) {

        // this is a substitution, create variant
        results[i].score = frequency_table[i][genome_pos].score;
        results[i].depth = frequency_table[i][genome_pos].score;
        results[i].ref[0] = nucleotide[ref_genome->data[genome_pos]];
        results[i].ref[1] = '\0';
        results[i].alt[0] = nucleotide[i];
        results[i].alt[1] = '\0';
    } else {
        results[i].score = 0;
        results[i].depth = 0;
    }
  }
  //printf("get_most_frequent_variant: genome_pos %lu, nucleotide max freq %d %f %c\n", genome_pos, nucId, max, nucId >= 0 ? nucleotide[nucId] : '-');

  // return results;
}

static void add_codependence_to_freq_table(struct frequency_info** frequency_table, struct variants_codependence_info_list* codependence_list, uint64_t position, uint8_t letter, float freq_update) {
    for (;codependence_list != NULL; codependence_list=codependence_list->next_list) {
        for (int i = 0; i<COD_LIST_SIZE && codependence_list->content[i].key != 0; i++) {
            uint8_t current_letter = codependence_list->content[i].key & 0x3;
            if (current_letter == letter) {
                uint8_t other_letter = (codependence_list->content[i].key >> 2) & 0x3;
                int64_t position_delta = codependence_list->content[i].key >> 4;
                frequency_table[other_letter][position + position_delta].freq += (codependence_list->content[i].codependence_count) * freq_update;
            }
        }
    }
}

#define POSITIVE_COD_INFLUENCE 0.0212905
#define NEGATIVE_COD_INFLUENCE -0.580356

//TODO here read frequency table and write vcf (take max of frequency table to find substitution if any)
void create_vcf()
{
    double start_time = my_clock();
    printf("%s:\n", __func__);

    FILE *vcf_file;
    int nb_variant = 0;
    char filename[1024];
    genome_t *ref_genome = genome_get();

    sprintf(filename, "%s_upvc.vcf", get_input_path());
    vcf_file = fopen(filename, "w");
    CHECK_FILE(vcf_file, filename);

    /* ####### START OF HEADER ####### */

    /* print vcf version (required) */
    fprintf(vcf_file, "##fileformat=VCFv4.3\n");

    /* print source of VCF file (this program) */
    fprintf(vcf_file, "##source=UPVC %s\n", VERSION);

    /* get the file date */
    char filedate[10];
    time_t mytime = time(NULL);
    strftime(filedate, 100, "%Y%d%m", localtime(&mytime));
    fprintf(vcf_file, "##fileDate=%s\n", filedate);

    /* print reference genome file name */
    fprintf(vcf_file, "##reference=%s.fasta\n", get_input_path());

    /* print the column names (fields are tab-delimited in VCF) */
    fprintf(vcf_file, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");

    /* ####### END OF HEADER ####### */

    struct frequency_info **frequency_table = get_frequency_table();
    struct variants_codependence_info_list** codependence_table = get_codependence_table(NULL);// giving a NULL pointer should only work if the table is already allocated; which should be the case
    uint32_t nb_pos_multiple_var = 0;

    /**
     * Dump debugging information: frequency table for a given set of positions
     **/
#if 0
    dbg_file = fopen("freq_debug.txt", "w");
    sub_file = fopen("subst.txt", "r");
    assert(sub_file);
    unsigned seq = 0;
    uint64_t pos = 0;
    static char nucleotide[4] = { 'A', 'C', 'T', 'G' };
    fprintf(dbg_file, "# seq pos ref-nucleotide frequencies:A C T G\n");
    while (EOF != fscanf(sub_file, "%u %lu\n", &seq, &pos))
    {
      assert(seq > 0);
      seq--;
      assert(pos > 0);
      pos--;
      for (int inc = -2; inc <= 2; inc++) {
        if(inc > 0 && pos + inc >= ref_genome->len_seq[seq]) continue;
        if(inc < 0 && pos < abs(inc)) continue;
        uint64_t genome_pos = ref_genome->pt_seq[seq] + pos + inc;
        if(genome_pos > genome_get()->fasta_file_size) {
          printf("WARNING: wrong genome position %lu. seq %u pos %lu inc %d\n", genome_pos, seq, pos, inc);
          continue;
        }
        fprintf(dbg_file, "%u %lu %c ", seq + 1, pos + inc + 1, nucleotide[ref_genome->data[genome_pos]]);
        float total = 0.0f;
        for(int m = 0; m < 5; ++m) {
            total += frequency_table[m][genome_pos].freq;
        }
        for(int m = 0; m < 4; ++m) {
            fprintf(dbg_file, "(%f %u %f %u) ",
            frequency_table[m][genome_pos].freq, frequency_table[m][genome_pos].score,
            frequency_table[m][genome_pos].freq * 100.0 / total, depth_filter(frequency_table[m][genome_pos].freq * 100.0 / total));
        }
        fprintf(dbg_file, "\n");
      }
    }
    fclose(dbg_file);
    fclose(sub_file);
#endif

#if DUMP_FREQUENCY_TABLE
    printf("dumping frequency table...\n");
    dump_freq_table(ref_genome, frequency_table);
    printf("table done dumping; starting variant calling...\n");
#endif

    unsigned int uncovered_nucleotides = 0;
    unsigned int badly_covered_nucleotides = 0;
    unsigned int well_covered_nucleotides = 0;
    unsigned int overly_covered_nucleotides = 0;
    unsigned int max_coverage = 0;
    unsigned int position_most_coverage = 0;
    unsigned int chromosome_most_coverage = 99999;
    uint64_t total_coverage = 0;
    uint64_t total_cov_squared = 0;
    variant_t variants_to_call[5];
    // First pass on the frequency table to take into account variant codependence
    LOG_INFO("doing first pass of vc\n");
    /* for each sequence in the genome */
    for (uint32_t seq_number = 0; seq_number < ref_genome->nb_seq; seq_number++) {
        /* for each position in the sequence */
        LOG_INFO("sequence %u\n", seq_number);
        for (uint64_t seq_position = 0; seq_position < ref_genome->len_seq[seq_number]; seq_position++) {
            get_most_frequent_variant(ref_genome, frequency_table, seq_number, seq_position, variants_to_call);
            for (uint8_t i=0; i<4; i++) {
                uint64_t genome_pos = ref_genome->pt_seq[seq_number] + seq_position;
                if (variants_to_call[i].depth) {
                    add_codependence_to_freq_table(frequency_table, codependence_table[genome_pos], genome_pos, i, POSITIVE_COD_INFLUENCE);
                } else if (frequency_table[i][genome_pos].score>0) {
                    add_codependence_to_freq_table(frequency_table, codependence_table[genome_pos], genome_pos, i, NEGATIVE_COD_INFLUENCE);
                }
            }
        }
    }

    free_codependence_chunks();

    LOG_INFO("doing second and final pass of vc\n");
    /* for each sequence in the genome */
    for (uint32_t seq_number = 0; seq_number < ref_genome->nb_seq; seq_number++) {
        /* for each position in the sequence */
        LOG_INFO("sequence %u\n", seq_number);
        for (uint64_t seq_position = 0; seq_position < ref_genome->len_seq[seq_number]; seq_position++) {
            
            get_most_frequent_variant(ref_genome, frequency_table, seq_number, seq_position, variants_to_call);
            unsigned int total_score = 0;
            total_score += frequency_table[0][ref_genome->pt_seq[seq_number] + seq_position].score;
            total_score += frequency_table[1][ref_genome->pt_seq[seq_number] + seq_position].score;
            total_score += frequency_table[2][ref_genome->pt_seq[seq_number] + seq_position].score;
            total_score += frequency_table[3][ref_genome->pt_seq[seq_number] + seq_position].score;
            total_coverage += total_score;
            total_cov_squared += total_score*total_score;
            //total_score += frequency_table[4][ref_genome->pt_seq[seq_number] + seq_position].score;
            if (total_score == 0) {
                    uncovered_nucleotides++;
            } else if (total_score < 10) {
                    badly_covered_nucleotides++;
            } else if (total_score < 90) {
                    well_covered_nucleotides++;
            } else {
                    overly_covered_nucleotides++;
                    if (total_score > max_coverage) {
                            max_coverage = total_score;
                            chromosome_most_coverage = seq_number;
                            position_most_coverage = seq_position;
                    }
            }
            int nb_var = 0;
            for(int i = 0; i < 5; ++i) {
              variant_t * var = &variants_to_call[i];
              if(var->depth) {
                // LOG_DEBUG("calling variant %d at %u:%lu, freq:%f\n", i, seq_number, seq_position, frequency_table[i][ref_genome->pt_seq[seq_number] + seq_position].freq);
                nb_variant += print_variant_tree(var, seq_number, seq_position, ref_genome, vcf_file) ? 1 : 0;
                nb_var++;
              }
              if(nb_var > 1)
                nb_pos_multiple_var++;
            }
        }
    }

    free_frequency_table();
    fclose(vcf_file);

    unsigned long total_nucleotides = overly_covered_nucleotides+well_covered_nucleotides+badly_covered_nucleotides+uncovered_nucleotides;
    printf("\tuncovered nucleotides: %lu (%lu.%lu%%)\n",
		    (long)uncovered_nucleotides,
		    (long)uncovered_nucleotides*100/total_nucleotides,
		    (long)uncovered_nucleotides*10000/total_nucleotides%100);
    printf("\tbadly covered nucleotides (less than 10 reads): %lu (%lu.%lu%%)\n",
		    (long)badly_covered_nucleotides,
		    (long)badly_covered_nucleotides*100/total_nucleotides,
		    (long)badly_covered_nucleotides*10000/total_nucleotides%100);
    printf("\twell covered nucleotides (10 to 90 reads): %lu (%lu.%lu%%)\n",
		    (long)well_covered_nucleotides,
		    (long)well_covered_nucleotides*100/total_nucleotides,
		    (long)well_covered_nucleotides*10000/total_nucleotides%100);
    printf("\toverly covered nucleotides (more than 90 reads): %lu (%lu.%lu%%)\n",
		    (long)overly_covered_nucleotides,
		    (long)overly_covered_nucleotides*100/total_nucleotides,
		    (long)overly_covered_nucleotides*10000/total_nucleotides%100);
    printf("\tmax coverage: %u reads\n", max_coverage);
    printf("\tmax coverage position: chr%u:%u\n", chromosome_most_coverage, position_most_coverage);
    printf("\ttotal coverage: %lu (eq %lu reads; or %lux coverage)\n", total_coverage, (long)total_coverage/SIZE_READ, (long)total_coverage/total_nucleotides);
   double mean = ((double)total_coverage) / (double) total_nucleotides;
    printf("\tmean cov: %f (std dev: %f)\n", mean, sqrt((double)total_cov_squared/(double)total_nucleotides - mean*mean));
    printf("\tnumber of variants: %d (multiple %d)\n", nb_variant, nb_pos_multiple_var);
    printf("\ttime: %lf s\n", my_clock() - start_time);
    fflush(stdout);
}
