#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <stdint.h>

#define NB_SEED ( 16*1024*1024 )

typedef struct index_seed {
        uint32_t nb_nbr;
        uint32_t offset;
        uint32_t num_dpu;
        struct index_seed *next;
} index_seed_t;

void save_index_seeds(index_seed_t **index_seed, char *seed_file)
{
        index_seed_t *seed;
        FILE *f = fopen(seed_file, "w");
        assert(f != NULL);

        for (unsigned int i = 0; i < NB_SEED; i++) {
                seed = index_seed[i];
                while (seed != NULL) {
                        fwrite(&i, sizeof(uint32_t), 1, f);
                        fwrite(seed, sizeof(uint32_t), 3, f);
                        seed = seed->next;
                }
        }

        fclose(f);
}

index_seed_t **load_index_seeds(char *seed_file)
{
        FILE *f = fopen(seed_file, "r");
        index_seed_t **index_seed;

        assert(f != NULL);

        printf("Loading index seeds\n");

        index_seed = (index_seed_t **) calloc(NB_SEED, sizeof(index_seed_t *));

        { /* First line is just a comment, skip */
                char line[512];
                assert(fgets(line, sizeof(line), f) != NULL);
        }

        {
                unsigned int seed_id = 0, dpu = 0, offset = 0, nb_nbr = 0;
                while (fscanf(f, "%u %u %u %u", &seed_id, &dpu, &offset, &nb_nbr) == 4) {
                        index_seed_t *seed = index_seed[seed_id];
                        if (seed == NULL) {
                                index_seed[seed_id] = seed = (index_seed_t *) malloc(sizeof(index_seed_t));
                        } else {
                                while (seed->next != NULL) {
                                        seed = seed->next;
                                }
                                seed->next = (index_seed_t *) malloc(sizeof(index_seed_t));
                                seed = seed->next;
                        }

                        seed->num_dpu = dpu;
                        seed->nb_nbr = nb_nbr;
                        seed->offset = offset;
                        seed->next = NULL;
                }
        }

        fclose(f);

        return index_seed;
}


int main(__attribute((unused)) int argc, char **argv)
{
        index_seed_t **index_seed = load_index_seeds(argv[1]);
        save_index_seeds(index_seed, argv[2]);
        return 0;
}
