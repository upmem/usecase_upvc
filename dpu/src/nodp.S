//
// Copyright 2016-2018 - Dominique Lavenier & UPMEM
//

#define COST_SUB    10

#define S1_ARG r0
#define S2_ARG r1
#define MAX_SCORE r2
#define SIZE_NEIGHBOUR_IN_BYTES r3

#define SCORE r0
#define S2 r1
#define S1 r4
#define S1_WORD r5
#define S2_WORD r6
// r7 and r8 unused for no reason other than historical
#define I r9
#define SIZE_NEIGHBOUR_IN_BYTES_MINUS_4 r10
#define S_XOR r11
#define TMP1 r12
#define S_TRANSLATED r12
#define TMP2 r13
#define V1 r14
#define V2 r15

    .macro CMP V1:req V2:req shift:req
    lsr TMP1, \V1, \shift
    and TMP2, \V2, (0xffffffff >> \shift)
    jeq TMP1, TMP2, .Lreturn_UINT_MAX
    .endm

    .macro CMP_PAIR V1:req V2:req shift:req
    CMP \V1, \V2, \shift
    CMP \V2, \V1, \shift
    .endm

    .text

    .section    .text.nodp,"ax",@progbits
    .global     nodp
    .global translation_table
    .type nodp,@function


nodp:
    sd r22, 0, d14
    move S1, S1_ARG
    move SCORE, 0
    move I, 0
.Lfor_i:
    // (9 to 29 instructions per cycle depending on xor result)
    lw S1_WORD, S1, 0
    lw S2_WORD, S2, 0
    add S1, S1, 4
    add S2, S2, 4
    xor S_XOR, S1_WORD, S2_WORD, z, .Lloopskip // skip k loop if xor is 0 as score doesn't need to be updated

.Lfor_k:
    // only run when S_XOR!=0 (8 instructions per cycle)
    and S_TRANSLATED, S_XOR, 0xff
    lbu S_TRANSLATED, S_TRANSLATED, translation_table
    add SCORE, SCORE, S_TRANSLATED

    jgtu SCORE, MAX_SCORE, .Lcheckoffset

    // loop
.Lloop_conditions:
    add I, I, 1
    jgeu I, SIZE_NEIGHBOUR_IN_BYTES, .Lend
    and zero, I, 0x3, z, .Lfor_i           // if (I%4==0) jump .Lfor_i
    lsr S_XOR, S_XOR, 8, nz, .Lfor_k       // S_XOR >>= 8; if (S_XOR!=0) jump .Lfor_k; //continue updating score only if necessary
                                           // otherwise, run loopskip

.Lloopskip:
    // k loop but without updating score (4 instructions total)
    and I, I, 0xfffffffc
    add I, I, 4
    jgeu I, SIZE_NEIGHBOUR_IN_BYTES, .Lend
    jump .Lfor_i


.Lcheckoffset:
    // first find the end
    and I, I, 0xfffffffc
    sub TMP1, SIZE_NEIGHBOUR_IN_BYTES_MINUS_4, I
    add S1, S1, TMP1
    add S2, S2, TMP1
    // Make sure an entire word can be loaded (go back into the read a little if needed)
    and S1, S1, 0xfffffffc
    and S2, S2, 0xfffffffc
    // load the last words
    lw S1_WORD, S1, 0
    lw S2_WORD, S2, 0
    CMP_PAIR S1_WORD, S2_WORD, 2
    CMP_PAIR S1_WORD, S2_WORD, 4
    CMP_PAIR S1_WORD, S2_WORD, 6
    CMP_PAIR S1_WORD, S2_WORD, 8
    CMP_PAIR S1_WORD, S2_WORD, 10
    jump .Lend

.Lreturn_UINT_MAX:
    move SCORE, 0xffffffff

.Lend:
    // TODO: process the last nucleotides when number of nucleotides is not a multiple of 4
    
    ld d14, r22, 0
    jump r23

    .size nodp, .Lend-nodp
    .section .stack_sizes,"o",@progbits,.text.nodp,unique,0
    .long nodp
    .byte 8
    .section    .text.nodp,"ax",@progbits
