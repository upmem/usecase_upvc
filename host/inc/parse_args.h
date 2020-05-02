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
 * @brief Get the goal of the run of the application.
 */
goal_t get_goal();

/**
 * @brief Get the number of DPUs to be use to compute.
 */
unsigned int get_nb_dpu();

/**
 * @brief Get the execution mode.
 */
bool get_simulation_mode();

bool get_no_filter();

/**
 * @brief Parse and validate the argument of the application.
 */
void validate_args(int argc, char **argv);

/**
 * @brief Free all the memory allocated during "validate_args".
 */
void free_args();

#endif /* __PARSE_ARGS_H__ */
