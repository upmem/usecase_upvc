/**
 * @Copyright (c) 2016-2018 - Dominique Lavenier & UPMEM
 */

#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <stdio.h>
#include <getopt.h>
#include <unistd.h>
#include "upvc.h"
#include "parse_args.h"

static char *prog_name = NULL;
static char *input_path = NULL;
static char *input_fasta = NULL;
static char *input_pe1 = NULL;
static char *input_pe2 = NULL;
static char *dpu_binary = NULL;
static bool simulation_mode = false;
static target_type_t target_type = target_type_unknown;
static goal_t goal = goal_unknown;
static nb_dpu_t nb_dpu = nb_dpu_unknown;

/**************************************************************************************/
/**************************************************************************************/
static void usage()
{
        ERROR_EXIT(24, "\nusage: %s -i <input_prefix> -d <number_of_dpus> [-s | -t <type> -g <goal> -b <dpu_binary>]\n"
                   "options:\n"
                   "\t-i\tInput prefix that will be used to find the inputs files\n"
                   "\t-d\tNumber of DPUs to use - value=128|256\n"
                   "\t-s\tSimulation mode (not compatible with -t -g and -b)\n"
                   "\t-t\tTarget type - values=hsim|fpga\n"
                   "\t-g\tGoal of the run - values=index|check|map\n"
                   "\t-b\tDPU binary to use"
                   ,
                   prog_name);
}

static void check_args()
{
        if (prog_name == NULL
            || input_path  == NULL
            || input_fasta == NULL
            || input_pe1   == NULL
            || input_pe2   == NULL
            || nb_dpu      == nb_dpu_unknown) {
                ERROR("missing option");
                usage();
        } else if (simulation_mode
                   && (target_type   != target_type_unknown
                       || dpu_binary != NULL)) {
                ERROR("simulation mode is not compatible with -t -g and -b option");
                usage();
        } else if (!simulation_mode
                   && (dpu_binary     == NULL
                       || target_type == target_type_unknown
                       || goal        == goal_unknown)) {
                ERROR("missing option");
                usage();
        }
}

/**************************************************************************************/
/**************************************************************************************/
static void verify_that_file_exists(const char *path)
{
        if (access(path, R_OK)) {
                ERROR_EXIT(25, "input file %s does not exist or is not readable (errno : %i)\n", path, errno);
        }
}

static char *alloc_input_file_name(const char *input_prefix, const char *input_suffix)
{
        char *input_file_name = (char *) calloc(strlen(input_prefix) + strlen(input_suffix) + 1, sizeof(char));
        sprintf(input_file_name, "%s%s", input_prefix, input_suffix);
        verify_that_file_exists(input_file_name);
        return input_file_name;
}

static void validate_inputs(const char *input_prefix)
{
        if (!(input_path     == NULL
              && input_fasta == NULL
              && input_pe1   == NULL
              && input_pe2   == NULL)) {
                ERROR("input option has been entered more than once");
                usage();
        } else {
                input_path = strdup(input_prefix);
                input_fasta = alloc_input_file_name(input_prefix, ".fasta");
                input_pe1   = alloc_input_file_name(input_prefix, "_PE1.fastq");
                input_pe2   = alloc_input_file_name(input_prefix, "_PE2.fastq");
        }
}

char *get_input_path()  { return input_path; }
char *get_input_fasta() { return input_fasta;}
char *get_input_pe1()   { return input_pe1;  }
char *get_input_pe2()   { return input_pe2;  }

/**************************************************************************************/
/**************************************************************************************/
static void validate_dpu_binary(const char *dpu_binary_str)
{
        if (dpu_binary != NULL) {
                ERROR("dpu binary option has been entered more than once");
                usage();
        } else {
                dpu_binary = strdup(dpu_binary_str);
                verify_that_file_exists(dpu_binary);
        }
}

char *get_dpu_binary() { return dpu_binary;}

/**************************************************************************************/
/**************************************************************************************/
static void validate_target_type(const char *target_type_str)
{
        if (target_type != target_type_unknown) {
                ERROR("target type option has been entered more than once");
                usage();
        } else if (strcmp(target_type_str, "fpga") == 0) {
                target_type = target_type_fpga;
        } else if (strcmp(target_type_str, "hsim") == 0) {
                target_type = target_type_hsim;
        } else {
                ERROR("unknown target type value");
                usage();
        }
}

target_type_t get_target_type() { return target_type;}

/**************************************************************************************/
/**************************************************************************************/
static void validate_goal(const char *goal_str)
{
        if (goal != goal_unknown) {
                ERROR("goal option has been entered more than once");
                usage();
        } else if (simulation_mode) {
                ERROR("goal is not compatible with simulation mode");
                usage();
        } else if (strcmp(goal_str, "index") == 0) {
                goal = goal_index;
        } else if (strcmp(goal_str, "check") == 0) {
                goal = goal_check;
        } else if (strcmp(goal_str, "map") == 0) {
                goal = goal_map;
        } else {
                ERROR("unknown goal value");
                usage();
        }
}

goal_t get_goal() { return goal;}

/**************************************************************************************/
/**************************************************************************************/
static void validate_nb_dpu(const char *nb_dpu_str)
{
        if (nb_dpu != nb_dpu_unknown) {
                ERROR("number of DPUs option has been entered more than once");
                usage();
        } else if (strcmp(nb_dpu_str, "128") == 0) {
                nb_dpu = nb_dpu_128;
        } else if (strcmp(nb_dpu_str, "256") == 0) {
                nb_dpu = nb_dpu_256;
        } else {
                ERROR("wrong number of DPUs");
                usage();
        }
}

nb_dpu_t get_nb_dpu() { return nb_dpu;}

/**************************************************************************************/
/**************************************************************************************/
static void validate_simulation_mode()
{
        simulation_mode = true;
        if (goal != goal_unknown) {
                ERROR("goal is not compatible with simulation mode");
                usage();
        }
        goal = goal_map;
}

bool get_simulation_mode() { return simulation_mode;}

/**************************************************************************************/
/**************************************************************************************/
void validate_args(int argc, char **argv)
{
        int opt;
        extern char *optarg;

        prog_name = strdup(argv[0]);
        while ((opt = getopt(argc, argv, "si:t:g:d:b:")) != -1) {
                switch (opt) {
                case 'i':
                        validate_inputs(optarg);
                        break;
                case 't':
                        validate_target_type(optarg);
                        break;
                case 'g':
                        validate_goal(optarg);
                        break;
                case 'd':
                        validate_nb_dpu(optarg);
                        break;
                case 'b':
                        validate_dpu_binary(optarg);
                        break;
                case 's':
                        validate_simulation_mode();
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
        free(input_fasta);
        free(input_pe1);
        free(input_pe2);
        free(dpu_binary);
}
