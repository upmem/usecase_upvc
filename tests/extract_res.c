#include <stdbool.h>
#include <stdio.h>
#include <string.h>

#define STRING_MAX_SIZE (1024)
#define FORMAT_R ("R:")

static bool read_log(FILE *fp, char *str)
{
    do {
        if (fgets(str, STRING_MAX_SIZE, fp) == NULL) {
            return false;
        }
    } while (strncmp(str, FORMAT_R, strlen(FORMAT_R)));

    return true;
}

static bool read_write_dpu(FILE *fp, char *str, unsigned int dpu_nb, char *file_name)
{
    char dpu_file_name[STRING_MAX_SIZE];
    unsigned int nb_res = 0;

    sprintf(dpu_file_name, "%s.%u", file_name, dpu_nb);
    printf("writing %s\n", dpu_file_name);

    FILE *write_fp = fopen(dpu_file_name, "w");

    do {
        nb_res++;
        fwrite(str, sizeof(char), strlen(str), write_fp);
        if (fgets(str, STRING_MAX_SIZE, fp) == NULL) {
            fclose(write_fp);
            return false;
        }
    } while (!strncmp(str, FORMAT_R, strlen(FORMAT_R)));

    printf("%u res for dpu %u\n", nb_res, dpu_nb);

    fclose(write_fp);
    return true;
}

int main(__attribute__((unused)) int argc, char **argv)
{
    unsigned int dpu_nb = 0;
    char str[STRING_MAX_SIZE];
    FILE *fp = fopen(argv[1], "r");

    while (true) {
        if (!read_log(fp, str)) {
            break;
        }
        if (!read_write_dpu(fp, str, dpu_nb, argv[1])) {
            break;
        }
        dpu_nb++;
    }
    printf("%u dpu read\n", dpu_nb);

    fclose(fp);
    return 0;
}
