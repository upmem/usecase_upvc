/**
 * @Copyright (c) 2016-2018 - Dominique Lavenier & UPMEM
 */

#ifndef __INTEGRATION_REQUEST_H__
#define __INTEGRATION_REQUEST_H__

/**
 * @param len  The size of a sequence, in bytes.
 * @return The corresponding number of symbols.
 */
static inline unsigned int nb_bytes_to_syms(unsigned int len) {
        return len << 2;
}

/**
 * @brief A request, as provided by the host.
 *
 * See dispatch.h. The structure definition does not include the very first bytes of neighbour.
 */
typedef struct {
        uint32_t offset;
        uint32_t count;
        uint32_t num;
        uint32_t nbr;
} request_t;

#define SIZEOF_REQUEST_HEADER (sizeof(request_t) - sizeof(uint32_t))

#endif /* __INTEGRATION_REQUEST_H__ */
