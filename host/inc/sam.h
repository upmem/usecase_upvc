/**
 * Copyright 2021 - A Moisson-Franckhauser & UPMEM
 */
#ifndef __SAM_H__
#define __SAM_H__

void open_sam_file();
//TODO : either reuse this code or delete it
//void write_sam_read(uint64_t genome_pos, uint8_t *code, int8_t *read);
void write_read_mapping(uint64_t genome_pos, uint8_t *code);
void close_sam_file();

#endif /* __SAM_H__ */
