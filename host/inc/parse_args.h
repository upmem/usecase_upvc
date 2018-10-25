#ifndef __PARSE_ARGS_H__
#define __PARSE_ARGS_H__

#include <stdbool.h>

typedef enum {target_type_unknown, target_type_hsim, target_type_fpga} target_type_t;
typedef enum{goal_unknown, goal_index, goal_check, goal_map} goal_t;
typedef enum{nb_dpu_unknown, nb_dpu_128 = 128, nb_dpu_256 = 256} nb_dpu_t;

char *get_input_path();
char *get_input_fasta();
char *get_input_pe1();
char *get_input_pe2();
char *get_dpu_binary();
target_type_t get_target_type();
goal_t get_goal();
nb_dpu_t get_nb_dpu();
bool get_simulation_mode();

void validate_args(int argc, char **argv);
void free_args();

#endif /* __PARSE_ARGS_H__ */
