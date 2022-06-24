#include <assert.h>
#include <stdio.h>
#include <string.h>

#define MAX_SIZE_LINE 512

int main(__attribute__((unused)) int argc, char **argv)
{
    FILE *chrall_file = fopen(argv[1], "r");
    unsigned int current_chr = 0;
    FILE *current_file = NULL;
    char str[MAX_SIZE_LINE];

    while (fgets(str, MAX_SIZE_LINE, chrall_file) != NULL) {
        if (str[0] == '>') {
            if (current_file != NULL)
                fclose(current_file);
            current_chr++;
            char filename[MAX_SIZE_LINE];
            sprintf(filename, "chr%u.fasta", current_chr);
            current_file = fopen(filename, "w");
            assert(current_file != NULL);
        }
        fwrite(str, strlen(str), 1, current_file);
    }

    fclose(current_file);
    fclose(chrall_file);

    return 0;
}
