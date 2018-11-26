/**
 * @Copyright (c) 2016-2018 - Dominique Lavenier & UPMEM
 */

#ifndef __PARSE_ARGS_H__
#define __PARSE_ARGS_H__

#include <stdbool.h>

typedef enum {target_type_unknown, target_type_hsim, target_type_fpga} target_type_t;
typedef enum{goal_unknown, goal_index, goal_check, goal_map} goal_t;
typedef enum{nb_dpu_unknown, nb_dpu_128 = 128, nb_dpu_256 = 256} nb_dpu_t;

/**
 * @brief Get the path where to store temporary and final file
 * (which is the path where the fasta file has been found).
 */
char *get_input_path();

/**
 * @brief Get the path to the fasta file.
 */
char *get_input_fasta();

/**
 * @brief Get the path to the first PE file.
 */
char *get_input_pe1();

/**
 * @brief Get the path to the second PE file.
 */
char *get_input_pe2();

/**
 * @brief Get the path to the binary to be load on the DPUs.
 */
char *get_dpu_binary();

/**
 * @brief Get the target type on which to compute.
 */
target_type_t get_target_type();

/**
 * @brief Get the goal of the run of the application.
 */
goal_t get_goal();

/**
 * @brief Get the number of DPUs to be use to compute.
 */
nb_dpu_t get_nb_dpu();

/**
 * @brief Get the execution mode.
 */
bool get_simulation_mode();

/**
 * @brief Get the number of DPUs to be use simultaneously per run.
 */
unsigned int get_nb_dpus_per_run();

/**
 * @brief Parse and validate the argument of the application.
 */
void validate_args(int argc, char **argv);

/**
 * @brief Free all the memory allocated during "validate_args".
 */
void free_args();

#endif /* __PARSE_ARGS_H__ */
