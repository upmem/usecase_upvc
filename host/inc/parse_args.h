/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#ifndef __PARSE_ARGS_H__
#define __PARSE_ARGS_H__

#include <stdbool.h>

typedef enum { goal_unknown, goal_index, goal_map } goal_t;

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
 * @brief Get the goal of the run of the application.
 */
goal_t get_goal();

/**
 * @brief Get the number of DPUs to be use to compute.
 */
unsigned int get_nb_dpu();

/**
 * @brief Set the number of DPUs to be use to compute.
 */
void set_nb_dpu(unsigned int val);

/**
 * @brief Get the execution mode.
 */
bool get_simulation_mode();

/**
 * @brief Get the number of DPUs to be use simultaneously per run.
 */
unsigned int get_nb_dpus_per_run();

/**
 * @brief Set the number of DPUs to be use simultaneously per run.
 */
void set_nb_dpus_per_run(unsigned int val);

/**
 * @brief Parse and validate the argument of the application.
 */
void validate_args(int argc, char **argv);

/**
 * @brief Free all the memory allocated during "validate_args".
 */
void free_args();

#endif /* __PARSE_ARGS_H__ */
