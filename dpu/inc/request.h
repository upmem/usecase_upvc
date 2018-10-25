/*
 * Copyright (c) 2014-2018 - uPmem
 */

#ifndef INTEGRATION_REQUEST_H
#define INTEGRATION_REQUEST_H

/**
 * @param len the size of a sequence, in bytes
 * @return The corresponding number of symbols.
 */
static inline unsigned int nr_bytes_to_syms(unsigned int len) {
    return len << 2;
}

/**
 * @brief A request, as provided by the host
 *
 * See dispatch.h. The structure definition does not include the very first bytes of neighbor.
 */
typedef struct {
    uint32_t offset;
    uint32_t count;
    uint32_t num;
    uint32_t nbr;
} request_t;

#define SIZEOF_REQUEST_HEADER (sizeof(request_t) - sizeof(uint32_t))
#endif // INTEGRATION_REQUEST_H
