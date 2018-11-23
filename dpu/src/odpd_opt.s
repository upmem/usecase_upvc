//
// Copyright (c) 2016-2018 - Dominique Lavenier & UPMEM
//

#include <sysdefs.s>

#define NB_DIAGS        15
#define COST_SUB        10
#define COST_GAPO       11
#define COST_GAPE       1
#define COST_INIT       99

#define LINE_SIZE_VAL   (6*4)
#define LINE_SIZE       ${LINE_SIZE_VAL}

#define d0off           ${  0*4 }
#define p0off           ${  2*4 }
#define p0off1          ${  3*4 }
#define q0off           ${  4*4 }
#define q0off1          ${  5*4 }
#define d1off           ${  6*4 }
#define p1off           ${  8*4 }
#define q1off           ${ 10*4 }

.text

.section    .text.odpd,"ax",@progbits
.global     $odpd
.global     $__M

// odpd is called from a C program, providing the arguments as follows:
//  - r0 = uint8_t *s1
//  - r1 = uint8_t *s2
//  - r2 = int max_score
//  - r3 = unsigned int nbr_sym_len
//  - r4 = tasklet id

odpd[                      \
    s1      = %r0         |\
    s2      = %r1         |\
    mxScore = %r2         |\
    len     = %r3         |\
    me      = %r4         |\
    matrix  = %r4         |\
    mpp     = %r5         |\
    mlp     = %r6         |\
    i       = %r7         |\
    j       = %r8         |\
    cost    = %r9         |\
    v99     = %r10        |\
    v1      = %r11        |\
    t0      = __R0        |\
    t1      = __R1        |\
    QP      = __R2        |\
    mnScore = __VAL_RET__ ]\
(
save_registers:
    sd      %r22    +0      %d0
    sd      %r22    +8      %d4
    sd      %r22    +16     %d6
    sd      %r22    +24     %d8
    sd      %r22    +32     %d10

get_tasklet_matrix:                                         // matrix = __M + 8 * tid * 3 * (len + 1)
    move    t0      +3                                      // Number of buffers (P, Q & D) = 3
    add     t1      len     +1                              // local matrix depth (len + 1)
    mul_ul_ul   t0      t0      me                          // 3 * tid
    mul_ul_ul   t0      t0      t1                          // 3 * tid * (len + 1)
    lsl     t0      t0      +3                              // 3 * tid * (len + 1) * 8
    lw      matrix  %zero   __M                             // load __M
    add     matrix  matrix  t0                              // matrix = __M + 3 * tid * (len + 1) * 8

init_matrix:
    move    cost    0                                       // Initialization value
    move    v99     COST_INIT                               // Initialization value
    move    mpp     matrix                                  // Matrix line (index j)
    move    t1      ${(NB_DIAGS/2+1)}                       // Loop counter
init_matrix_loop:
    sw      @mpp    +d0off  cost                            // Set D(0, j)
    sw      @mpp    +p0off  v99                             // Set P(0, j)
    sw      @mpp    +q0off  v99                             // Set Q(0, j)
    add     cost    cost    +COST_SUB                       //
    add     mpp     mpp     +LINE_SIZE                      //
    add     t1      t1      -1:8   ?nz @init_matrix_loop    //
init_matrix_tail:
    sw      @matrix +p0off1 v99                             // Set P(1, 0)
    sw      @matrix +q0off1 v99                             // Set Q(1, 0)

phase1:
    move    cost    COST_SUB                                // Initialization value
    move    i       1                                       // Initialization value
    add     mpp     matrix  +4                              // Matrix line (index pp,i)
    add     mlp     matrix  +0                              // Matrix line (index lp,i)
phase1_iloop:
    sw      @mpp    +d0off  cost                            // Set D(pp, 0)
    move    mnScore v99                                     // Initialization value

phase1_v1:
    add     t0      i       -1                              // Probably requires optimization here
    ror     t0      t0      2                               //
    add     t1      s1      t0                              //
    lbu     t1      @t1     +0                              //
    lsr     t0      t0      29                              //
    lsr     v1      t1      t0                              //

    move    j       0                                       // Initialization value
    add     i       i       +${NB_DIAGS/2}                  // j loop bound
    lw      t0      @mpp    +d0off                          // Get first D(pp, j-1)
phase1_jloop:
    add     t0      t0      +COST_GAPO                      //
    lw      t1      @mpp    +p0off                          // Get P(pp, j-1)
    add     t1      t1      +COST_GAPE                      //
    jltu    t0      t1      @phase1_min0                    // Min is D(pp, j-1)+COST_GAPO
    move    t0      t1                                      // Min is P(pp, j-1)+COST_GAPE
phase1_min0:
    sw      @mpp    +p1off  t0                              // Set P(pp, j)
    lw      QP      @mlp    +d1off                          // Get D(lp, j)
    add     QP      QP      +COST_GAPO                      //
    lw      t1      @mlp    +q1off                          // Get Q(lp, j)
    add     t1      t1      +COST_GAPE                      //
    jltu    QP      t1      @phase1_min1                    // Min is D(lp, j)+COST_GAPO
    move    QP      t1                                      // Min is Q(lp, j)+COST_GAPE
phase1_min1:
    sw      @mpp    +q1off  QP                              // Set Q(pp, j)
    jltu    QP      t0      @phase1_min2                    // Min is Q(pp, j)
    move    QP      t0                                      // Min is P(pp, j)
phase1_min2:

phase1_v2:
    ror     t0      j       2                               //
    add     t1      s2      t0                              //
    lbu     t1      @t1     +0                              //
    lsr     t0      t0      29                              //
    lsr     t1      t1      t0                              //
    xor     t1      t1      v1                              // base difference

phase1_dcost:
    lw      t0      @mlp    +d0off                          // Get D(lp, j-1)
    and     t1      t1      3       ?z @phase1_nodiff       //
    add     t0      t0      +COST_SUB                       //
phase1_nodiff:
    jltu    t0      QP      @phase1_min3                    //
    move    t0      QP                                      //
phase1_min3:
    sw      @mpp    +d1off  t0                              // Set D(pp, j), upcoming D(pp, j-1)
    jltu    mnScore t0      @phase1_min4                    // Min score saved
    move    mnScore t0                                      // Min score updated
phase1_min4:

phase1_jnext:
    add     mpp     mpp     +LINE_SIZE                      // Loop increment
    add     mlp     mlp     +LINE_SIZE                      // Loop increment
    add     j       j       +1                              // Loop increment
    jltu    j       i       @phase1_jloop                   // Next iteration

phase1_earlyexit:
    sw      @mpp    +q1off  v99                             // Set Q(pp, j)
    sw      @mpp    +d1off  v99                             // Set D(pp, j)
    jltu    mxScore mnScore @restore_registers              // Early exit

phase1_inext:
    add     cost    cost    +COST_SUB                       // Loop increment
    add     i       i       ${1-NB_DIAGS/2}                 // Loop increment
    and     mpp     i       1       ?z @phase1_odd          //
    add     mpp     matrix  +4                              // Switch vector
    add     mlp     matrix  +0      ?true @phase1_itest     // Switch vector
phase1_odd:
    add     mpp     matrix  +0                              // Switch vector
    add     mlp     matrix  +4                              // Switch vector
phase1_itest:
    sub     t0     i ${NB_DIAGS/2+1}   ?ltu @phase1_iloop   // Next iteration

phase2:
phase2_iloop:
    move    mnScore v99                                     // Initialization value
    sw      @mpp    +p0off  v99                             // Set P(pp, j)
    sw      @mpp    +d0off  v99                             // Set D(pp, j)

phase2_v1:
    add     t0      i       -1                              // Probably requires optimization here
    ror     t0      t0      2                               //
    add     t1      s1      t0                              //
    lbu     t1      @t1     +0                              //
    lsr     t0      t0      29                              //
    lsr     v1      t1      t0                              //

    add     j       i       ${-NB_DIAGS/2-1}                // Initialization value
    add     i       i       +${NB_DIAGS/2}                  // j loop bound
    lw      t0      @mpp    +d0off                          // Get first D(pp, j-1)
phase2_jloop:
    add     t0      t0      +COST_GAPO                      //
    lw      t1      @mpp    +p0off                          // Get P(pp, j-1)
    add     t1      t1      +COST_GAPE                      //
    jltu    t0      t1      @phase2_min0                    // Min is D(pp, j-1)+COST_GAPO
    move    t0      t1                                      // Min is P(pp, j-1)+COST_GAPE
phase2_min0:
    sw      @mpp    +p1off  t0                              // Set P(pp, j)
    lw      QP      @mlp    +d1off                          // Get D(lp, j)
    add     QP      QP      +COST_GAPO                      //
    lw      t1      @mlp    +q1off                          // Get Q(lp, j)
    add     t1      t1      +COST_GAPE                      //
    jltu    QP      t1      @phase2_min1                    // Min is D(lp, j)+COST_GAPO
    move    QP      t1                                      // Min is Q(lp, j)+COST_GAPE
phase2_min1:
    sw      @mpp    +q1off  QP                              // Set Q(pp, j)
    jltu    QP      t0      @phase2_min2                    // Min is Q(pp, j)
    move    QP      t0                                      // Min is P(pp, j)
phase2_min2:

phase2_v2:
    ror     t0      j       2                               //
    add     t1      s2      t0                              //
    lbu     t1      @t1     +0                              //
    lsr     t0      t0      29                              //
    lsr     t1      t1      t0                              //
    xor     t1      t1      v1                              // base difference

phase2_dcost:
    lw      t0      @mlp    +d0off                          // Get D(lp, j-1)
    and     t1      t1      3       ?z @phase2_nodiff       //
    add     t0      t0      +COST_SUB                       //
phase2_nodiff:
    jltu    t0      QP      @phase2_min3                    //
    move    t0      QP                                      //
phase2_min3:
    sw      @mpp    +d1off  t0                              // Set D(pp, j), upcoming D(pp, j-1)
    jltu    mnScore t0      @phase2_min4                    // Min score saved
    move    mnScore t0                                      // Min score updated
phase2_min4:

phase2_jnext:
    add     mpp     mpp     +LINE_SIZE                      // Loop increment
    add     mlp     mlp     +LINE_SIZE                      // Loop increment
    add     j       j       +1                              // Loop increment
    jltu    j       i       @phase2_jloop                   // Next iteration

phase2_earlyexit:
    sw      @mpp    +q1off  v99                             // Set Q(pp, j)
    sw      @mpp    +d1off  v99                             // Set D(pp, j)
    jltu    mxScore mnScore @restore_registers              // Early exit

phase2_inext:
    add     mpp     mpp     ${-LINE_SIZE_VAL*2*(NB_DIAGS/2)}           // Loop increment
    add     mlp     mlp     ${-LINE_SIZE_VAL*2*(NB_DIAGS/2)}           // Loop increment
    xor     mpp     mpp     4                               // Switch vector
    xor     mlp     mlp     4                               // Switch vector
    add     i       i       ${1-NB_DIAGS/2}                 // Loop increment
    add     t0      len     ${-(NB_DIAGS/2)}                // Loop bound (to be optimized)
    jltu    i       t0      @phase2_iloop                   // Next iteration

phase3:
    move    mnScore v99                                     // Initialization value
    add     matrix  mpp     +LINE_SIZE                      // Next start

phase3_iloop:
    sw      @mpp    +p0off  v99                             // Set P(pp, j)
    sw      @mpp    +d0off  v99                             // Set D(pp, j)

phase3_v1:
    add     t0      i       -1                              // Probably requires optimization here
    ror     t0      t0      2                               //
    add     t1      s1      t0                              //
    lbu     t1      @t1     +0                              //
    lsr     t0      t0      29                              //
    lsr     v1      t1      t0                              //

    add     j       i       ${-NB_DIAGS/2-1}                // Initialization value
    lw      t0      @mpp    +d0off                          // Get first D(pp, j-1)
phase3_jloop:
    add     t0      t0      +COST_GAPO                      //
    lw      t1      @mpp    +p0off                          // Get P(pp, j-1)
    add     t1      t1      +COST_GAPE                      //
    jltu    t0      t1      @phase3_min0                    // Min is D(pp, j-1)+COST_GAPO
    move    t0      t1                                      // Min is P(pp, j-1)+COST_GAPE
phase3_min0:
    sw      @mpp    +p1off  t0                              // Set P(pp, j)
    lw      QP      @mlp    +d1off                          // Get D(lp, j)
    add     QP      QP      +COST_GAPO                      //
    lw      t1      @mlp    +q1off                          // Get Q(lp, j)
    add     t1      t1      +COST_GAPE                      //
    jltu    QP      t1      @phase3_min1                    // Min is D(lp, j)+COST_GAPO
    move    QP      t1                                      // Min is Q(lp, j)+COST_GAPE
phase3_min1:
    sw      @mpp    +q1off  QP                              // Set Q(pp, j)
    jltu    QP      t0      @phase3_min2                    // Min is Q(pp, j)
    move    QP      t0                                      // Min is P(pp, j)
phase3_min2:

phase3_v2:
    ror     t0      j       2                               //
    add     t1      s2      t0                              //
    lbu     t1      @t1     +0                              //
    lsr     t0      t0      29                              //
    lsr     t1      t1      t0                              //
    xor     t1      t1      v1                              // base difference

phase3_dcost:
    lw      t0      @mlp    +d0off                          // Get D(lp, j-1)
    and     t1      t1      3       ?z @phase3_nodiff       //
    add     t0      t0      +COST_SUB                       //
phase3_nodiff:
    jltu    t0      QP      @phase3_min3                    //
    move    t0      QP                                      //
phase3_min3:
    sw      @mpp    +d1off  t0                              // Set D(pp, j), upcoming D(pp, j-1)

phase3_jnext:
    add     mpp     mpp     +LINE_SIZE                      // Loop increment
    add     mlp     mlp     +LINE_SIZE                      // Loop increment
    add     j       j       +1                              // Loop increment
    jltu    j       len     @phase3_jloop                   // Next iteration

phase3_min4:
    jltu    mnScore t0      @phase3_inext                   // Min is mnScore
    move    mnScore t0                                      // Min is D(pp, j-1)

phase3_inext:
    sub     t0      mpp     matrix                          // Offset to next line
    sub     mlp     mlp     t0                              // Loop increment
    xor     mpp     matrix  4                               // Switch vector
    xor     mlp     mlp     4                               // Switch vector
    add     matrix  mpp     +LINE_SIZE                      // Next start
    add     i       i       +1                              // Loop increment
    jleu    i       len     @phase3_iloop                   // Next iteration

phase4:
    xor     mpp     mpp     4                               // Switch vector
    add     j       len     ${-(NB_DIAGS/2)}                // Initialization value
phase4_jloop:
    lw      t0      @mpp    +d0off                          // Get D(pp,j)
    jltu    mnScore t0      @phase4_jnext                   // Min is mnScore
    move    mnScore t0                                      // Min is D(pp,j)
phase4_jnext:
    add     mpp     mpp     +LINE_SIZE                      // Loop increment
    add     j       j       +1                              // Loop increment
    jleu    j       len     @phase4_jloop                   // Next iteration

restore_registers:
    ld      %d0     %r22    +0
    ld      %d4     %r22    +8
    ld      %d6     %r22    +16
    ld      %d8     %r22    +24
    ld      %d10    %r22    +32
    return
) // odpd
