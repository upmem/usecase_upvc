/**
 * Copyright 2016-2019 - Dominique Lavenier & UPMEM
 */

#ifndef __VMI_H__
#define __VMI_H__

#include <stdbool.h>

typedef struct {
    int handle;
} vmi_t;

bool vmi_open(unsigned int dpu_id, vmi_t *vmi);

void vmi_close(vmi_t *vmi);

bool vmi_write(vmi_t *vmi, long offs, const void *buffer, size_t nb_bytes);

#endif /* __VMI_H__ */
