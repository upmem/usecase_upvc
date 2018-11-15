/**
 * @Copyright (c) 2016-2018 - Dominique Lavenier & UPMEM
 */

#ifndef __VMI_H__
#define __VMI_H__

#include <stdbool.h>

/**
 * @brief Virtual Memory Image: a swap system to save memory.
 *
 * One VMI is declared with a name, defining the file that will mirror memory contents.
 * It provides the following functions:
 *   - Append a buffer to the end of memory
 *   - Copy a buffer at a specific location in memory
 *   - Copy the memory contents from
 *
 */

/**
 * @brief A virtual memory image descriptor.
 */
typedef struct {
        char *file_name;
        size_t mem_size;
} vmi_t;

/**
 * @brief Creates a new VMI
 *
 * @param name  VMI assigned name, from which the file name is guessed.
 * @param vmi   the VMI structure, initialized by this function.
 *
 * @return False if the memory could not be created: cannot access file, or memory exhausted.
 */
bool vmi_create(const char *name, vmi_t *vmi);

/**
 * @brief Deletes a VMI
 *
 * The VMI file remains on disk.
 *
 * @param vmi The deleted structure.
 */
void vmi_delete(vmi_t *vmi);

/**
 * @param vmi  A virtual memory.
 *
 * @return The current memory size, in bytes.
 */
static inline size_t vmi_size(const vmi_t *vmi) {
        return vmi->mem_size;
}

/**
 * @brief Pushes a buffer at the end of a VMI
 *
 * @param vmi       The virtual memory.
 * @param buffer    The buffer.
 * @param nb_bytes  The buffer size, in bytes.
 *
 * @return False if the operation failed for some unexpected reason.
 */
bool vmi_push(vmi_t *vmi, const void *buffer, size_t nb_bytes);

/**
 * @brief Copies a buffer at an exact location in a VMI.
 *
 * @param vmi       The virtual memory.
 * @param offs      Target offset in memory (must be less than the current memory size).
 * @param buffer    The buffer.
 * @param nb_bytes  The buffer size, in bytes.
 *
 * @return False if the operation failed for some unexpected reason.
 */
bool vmi_write(vmi_t *vmi, long offs, const void *buffer, size_t nb_bytes);

/**
 * @brief Reads the contents of a VMI.
 *
 * @param vmi       The virtual memory.
 * @param buffer    Target buffer.
 * @param nb_bytes  The buffer size, in bytes.
 *
 * @return The number of bytes actually read.
 */
size_t vmi_read(vmi_t *vmi, void *buffer, size_t nb_bytes);

#endif /* __VMI_H__ */
