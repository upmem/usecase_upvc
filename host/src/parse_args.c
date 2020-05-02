/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#define _GNU_SOURCE
#define _POSIX_C_SOURCE 200809L
#include <assert.h>
#include <errno.h>
#include <getopt.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <dpu.h>

#include "parse_args.h"
#include "upvc.h"

static char *prog_name = NULL;
static char *input_path = NULL;
static bool simulation_mode = false;
static bool no_filter = false;
static goal_t goal = goal_unknown;
static unsigned int nb_dpu = DPU_ALLOCATE_ALL;

/**************************************************************************************/
/**************************************************************************************/
static void usage()
{
    ERROR_EXIT(ERR_USAGE,
        "\nusage: %s -i <input_prefix> -g <goal> [-s | -n <number_of_dpus>] \n"
        "options:\n"
        "\t-i\tInput prefix that will be used to find the inputs files\n"
        "\t-g\tGoal of the run - values=index|map\n"
        "\t-s\tSimulation mode (not compatible with -n)\n"
        "\t-n\tNumber of DPUs to use when not in simulation mode\n",
        prog_name);
}

static void check_args()
{
    if (simulation_mode && nb_dpu != DPU_ALLOCATE_ALL) {
        ERROR("-n is not compatible with simulation mode");
        usage();
    }
    if (prog_name == NULL || input_path == NULL || goal == goal_unknown) {
        ERROR("missing option");
        usage();
    }
    if (goal == goal_index) {
        if (nb_dpu == DPU_ALLOCATE_ALL) {
            ERROR("missing option (number of dpus)");
            usage();
        } else if (nb_dpu == 0) {
            ERROR("cannot index for 0 dpus");
            usage();
        }
    }
}

static void check_permission()
{
    if (access(".", R_OK | W_OK)) {
        ERROR_EXIT(ERR_CURRENT_FOLDER_PERMISSIONS, "'%s' does not have read write permission on the current folder\n", prog_name);
    }
}

/**************************************************************************************/
/**************************************************************************************/

static void validate_inputs(const char *input_prefix)
{
    if (input_path != NULL) {
        ERROR("input option has been entered more than once");
        usage();
    } else {
        input_path = strdup(input_prefix);
        assert(input_path != NULL);
    }
}

char *get_input_path() { return input_path; }

/**************************************************************************************/
/**************************************************************************************/
static void validate_goal(const char *goal_str)
{
    if (goal != goal_unknown) {
        ERROR("goal option has been entered more than once");
        usage();
    } else if (strcmp(goal_str, "index") == 0) {
        goal = goal_index;
    } else if (strcmp(goal_str, "map") == 0) {
        goal = goal_map;
    } else {
        ERROR("unknown goal value");
        usage();
    }
}

goal_t get_goal() { return goal; }

/**************************************************************************************/
/**************************************************************************************/
static void validate_nb_dpu(const char *nb_dpu_str)
{
    if (nb_dpu != DPU_ALLOCATE_ALL) {
        ERROR("number of DPUs per run option has been entered more than once");
        usage();
    }
    nb_dpu = (unsigned int)atoi(nb_dpu_str);
}

unsigned int get_nb_dpu() { return nb_dpu; }

/**************************************************************************************/
/**************************************************************************************/
static void validate_simulation_mode() { simulation_mode = true; }

bool get_simulation_mode() { return simulation_mode; }

/**************************************************************************************/
/**************************************************************************************/
static void validate_no_filter() { no_filter = true; }

bool get_no_filter() { return no_filter; }

/**************************************************************************************/
/**************************************************************************************/
void validate_args(int argc, char **argv)
{
    int opt;
    extern char *optarg;

    prog_name = strdup(argv[0]);
    check_permission();

    while ((opt = getopt(argc, argv, "fsi:g:n:")) != -1) {
        switch (opt) {
        case 'i':
            validate_inputs(optarg);
            break;
        case 'g':
            validate_goal(optarg);
            break;
        case 's':
            validate_simulation_mode();
            break;
        case 'n':
            validate_nb_dpu(optarg);
            break;
        case 'f':
            validate_no_filter();
            break;
        default:
            ERROR("unknown option");
            usage();
        }
    }
    check_args();
}

void free_args()
{
    free(prog_name);
    free(input_path);
}
