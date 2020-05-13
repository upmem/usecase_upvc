//
// Copyright 2016-2018 - Dominique Lavenier & UPMEM
//

#define COST_SUB    10
#define SIZE_NEIGHBOUR_IN_BYTES 26

#define S1_ARG r0
#define S2_ARG r1
#define DELTA r2
#define MAX_SCORE r3

#define SCORE r0
#define S2 r1
#define S1 r2
#define S1_ACC_H r4
#define S1_ACC_L r5
#define S2_ACC_H r6
#define S2_ACC_L r7
#define I r8
#define S_XOR r9
#define TMP1 r10
#define S_TRANSLATED r10
#define TMP2 r11
#define V1 r12
#define V2 r13

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
    move S1, S1_ARG
    move SCORE, 0
    move I, 0
    lw S1_ACC_H, S1, 0
    lw S2_ACC_H, S2, 0
.Lfor_i:
    move S1_ACC_L, S1_ACC_H
    move S2_ACC_L, S2_ACC_H
    add S1, S1, 4
    add S2, S2, 4
    lw S1_ACC_H, S1, 0
    lw S2_ACC_H, S2, 0
    xor S_XOR, S1_ACC_L, S2_ACC_L
.Lfor_k:
    and S_TRANSLATED, S_XOR, 0xff
    lbu S_TRANSLATED, S_TRANSLATED, translation_table
    add SCORE, SCORE, S_TRANSLATED
    jleu S_TRANSLATED, COST_SUB, .Lupdate_score

    jgeu I, (SIZE_NEIGHBOUR_IN_BYTES - 4), .Lupdate_score

    and TMP2, I, 0x3
    add TMP2, TMP2, 1
    lsl TMP2, TMP2, 3

    lsrx TMP1, S1_ACC_H, TMP2, sh32, .LS1_J32
    lsr V1, S1_ACC_L, TMP2
    add V1, V1, TMP1, true, .LS1
.LS1_J32:
    move V1, S1_ACC_H
.LS1:
    lsrx TMP1, S2_ACC_H, TMP2, sh32, .LS2_J32
    lsr V2, S2_ACC_L, TMP2
    add V2, V2, TMP1, true, .LS2
.LS2_J32:
    move V2, S2_ACC_H
.LS2:

    CMP_PAIR V1, V2, 2
    CMP_PAIR V1, V2, 4
    CMP_PAIR V1, V2, 6
    CMP_PAIR V1, V2, 8

.Lupdate_score:
    jgtu SCORE, MAX_SCORE, .Lend

    // loop
    add I, I, 1
    jgeu I, SIZE_NEIGHBOUR_IN_BYTES, .Lend
    and zero, I, 0x3, z, .Lfor_i
    lsr S_XOR, S_XOR, 8, true, .Lfor_k

.Lreturn_UINT_MAX:
    move SCORE, 0xffffffff
.Lend:
    jump r23

    .size nodp, .Lend-nodp
    .section .stack_sizes,"o",@progbits,.text.nodp,unique,0
    .long nodp
    .byte 8
    .section    .text.nodp,"ax",@progbits
