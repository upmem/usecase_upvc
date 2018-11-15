/**
 * @Copyright (c) 2016-2018 - Dominique Lavenier & UPMEM
 */

#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>

#include "upvc.h"
#include "vmi.h"

#define VMI_EXT (".vmi")

static bool create_or_clear_vmi_file(const char *file_name)
{
        FILE *file = fopen(file_name, "w");
        if (file == NULL) {
                ERROR("VMI: could not create file for writing");
                return false;
        }

        fclose(file);
        return true;
}

bool vmi_create(const char *name, vmi_t *vmi)
{
        vmi->file_name = (char *) malloc(strlen(name) + strlen(VMI_EXT) + 1);
        if (vmi->file_name == NULL) {
                ERROR("VMI: no memory left!");
                return false;
        }

        sprintf(vmi->file_name, "%s%s", name, VMI_EXT);
        vmi->mem_size = 0;

        if (!create_or_clear_vmi_file(vmi->file_name)) {
                vmi_delete(vmi);
                return false;
        }

        return true;
}

void vmi_delete(vmi_t *vmi)
{
        free(vmi->file_name);
}

bool vmi_push(vmi_t *vmi, const void *buffer, size_t nb_bytes)
{
        size_t written;
        FILE *file = fopen(vmi->file_name, "a");
        if (file == NULL) {
                ERROR("VMI: push failed: could not re-open file");
                return false;
        }

        written = fwrite(buffer, 1, nb_bytes, file);
        if (written != nb_bytes) {
                ERROR("VMI: push failed: wrote wrong number of bytes");
                return false;
        }

        fclose(file);
        vmi->mem_size += nb_bytes;
        return true;
}

bool vmi_write(vmi_t *vmi, long offs, const void *buffer, size_t nb_bytes)
{
        int handle = open(vmi->file_name, O_WRONLY | O_CREAT, 644);
        if (handle == -1) {
                ERROR("VMI: write failed: could not re-open file");
                return false;
        }

        lseek(handle, offs, SEEK_SET);
        write(handle, buffer, nb_bytes);
        close(handle);
        vmi->mem_size += nb_bytes;

        return true;
}

size_t vmi_read(vmi_t *vmi, void *buffer, size_t nb_bytes)
{
        size_t result;
        FILE *file = fopen(vmi->file_name, "r");
        if (file == NULL) {
                ERROR("VMI: failed to open the file for reading");
                return 0L;
        }

        result = fread(buffer, 1, nb_bytes, file);
        fclose(file);
        return result;
}
