/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#include <fcntl.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "mram_dpu.h"
#include "upvc.h"
#include "vmi.h"

bool vmi_open(unsigned int dpu_id, vmi_t *vmi)
{
    vmi->handle = open(make_mram_file_name(dpu_id), O_WRONLY | O_CREAT, S_IRUSR | S_IRGRP | S_IROTH);
    if (vmi->handle == -1) {
        ERROR("VMI: open failed");
        return false;
    }

    return true;
}

void vmi_close(vmi_t *vmi) { close(vmi->handle); }

bool vmi_write(vmi_t *vmi, long offs, const void *buffer, size_t nb_bytes)
{
    lseek(vmi->handle, offs, SEEK_SET);
    size_t write_size = write(vmi->handle, buffer, nb_bytes);
    if (write_size != nb_bytes) {
        ERROR("VMI: did not write the expected number of bytes");
        return false;
    }

    return true;
}
