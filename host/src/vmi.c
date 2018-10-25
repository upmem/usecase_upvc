#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>

#include "upvc.h"
#include "vmi.h"

static bool create_or_clear_vmi_file(const char *file_name) {
    FILE *file = fopen(file_name, "wb");
    if (file == NULL) {
        ERROR("VMI: could not create file for writing");
        return false;
    }

    fclose(file);
    return true;
}

bool vmi_create(const char *name, vmi_t *vmi) {
    vmi->file_name = (char *) malloc(strlen(name) + strlen(".vmi") + 1);
    if (vmi->file_name == NULL) {
        ERROR("VMI: no memory left!");
        return false;
    }

    (void) sprintf(vmi->file_name, "%s.vmi", name);
    vmi->mem_size = 0;

    if (!create_or_clear_vmi_file(vmi->file_name)) {
        vmi_delete(vmi);
        return false;
    }

    return true;
}

void vmi_delete(vmi_t *vmi) {
    free(vmi->file_name);
}

bool vmi_push(vmi_t *vmi, const void *buffer, size_t nr_bytes) {
    FILE *file = fopen(vmi->file_name, "ab");
    if (file == NULL) {
        ERROR("VMI: push failed: could not re-open file");
        return false;
    }

    size_t written = fwrite(buffer, 1, nr_bytes, file);
    if (written != nr_bytes) {
        ERROR("VMI: push failed: wrote wrong number of bytes");
        return false;
    }

    fclose(file);
    vmi->mem_size += nr_bytes;
    return true;
}

bool vmi_write(vmi_t *vmi, long offs, const void *buffer, size_t nr_bytes) {
    int handle;
    handle = open(vmi->file_name, O_WRONLY | O_CREAT, 644);
    if (handle == -1) {
        ERROR("VMI: write failed: could not re-open file");
        return false;
    }

    lseek(handle, offs, SEEK_SET);
    write(handle, buffer, nr_bytes);
    close(handle);
    vmi->mem_size += nr_bytes;

    return true;
}

size_t vmi_read(vmi_t *vmi, long offs, void *buffer, size_t nr_bytes) {
    FILE *file = fopen(vmi->file_name, "rb");
    if (file == NULL) {
        ERROR("VMI: failed to open the file for reading");
        return 0L;
    }

    fseek(file, offs, SEEK_SET);

    size_t result = fread(buffer, 1, nr_bytes, file);
    fclose(file);
    return result;
}
