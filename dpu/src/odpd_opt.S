//
// Copyright 2016-2018 - Dominique Lavenier & UPMEM
//

#define NB_DIAGS    15
#define COST_SUB    10
#define COST_GAPO   11
#define COST_GAPE   1
#define COST_INIT   999

#define LINE_SIZE   (  6*4 )
#define d0off       (  0*4 )
#define p0off       (  2*4 )
#define p0off1      (  3*4 )
#define q0off       (  4*4 )
#define q0off1      (  5*4 )
#define d1off       (  6*4 )
#define p1off       (  8*4 )
#define q1off       ( 10*4 )

#define mnScore     r0
#define s1          r4
#define s2          r1
#define mxScore     r2
#define len         r3
#define matrix      r5
#define mpp         r6
#define mlp         r7
#define i           r8
#define j           r9
#define cost        r10
#define v99         r11
#define v1          r12
#define QP          r13
#define t0          r14
#define t1          r15

.text

.section    .text.odpd,"ax",@progbits
.global     odpd
.global     __M
.type odpd,@function

// odpd is called from a C program, providing the arguments as follows:
//  - r0 = uint8_t *s1
//  - r1 = uint8_t *s2
//  - r2 = int max_score
//  - r3 = unsigned int nbr_sym_len

odpd:
    sd         r22,     0,       d14
    move       s1,      r0
.Lget_tasklet_matrix:                                                 // matrix = __M + 8 * tid * 3 * (len + 2)
    move       t0,      +3                                            // Number of buffers (P, Q & D) = 3
    add        t1,      len,     +2                                   // local matrix depth (len + 2)
    mul_ul_ul  t0,      id,      t0                                   // 3 * tid
    mul_ul_ul  t0,      t0,      t1                                   // 3 * tid * (len + 2)
    lsl        t0,      t0,      +3                                   // 3 * tid * (len + 2) * 8
    add        matrix,  t0,      __M                                  // matrix = __M + 3 * tid * (len + 2) * 8

.Linit_matrix:
    move       cost,    0                                             // Initialization value
    move       v99,     COST_INIT                                     // Initialization value
    move       mpp,     matrix                                        // Matrix line (index j)
    move       t1,      (NB_DIAGS/2+1)                                // Loop counter
.Linit_matrix_loop:
    sw         mpp,     +d0off,  cost                                 // Set D(0, j)
    sw         mpp,     +p0off,  v99                                  // Set P(0, j)
    sw         mpp,     +q0off,  v99                                  // Set Q(0, j)
    add        cost,    cost,    +COST_SUB                            //
    add        mpp,     mpp,     +LINE_SIZE                           //
    add        t1,      t1,      -1, ?pl, .Linit_matrix_loop          //
.Linit_matrix_tail:
    sw         matrix,  +p0off1, v99                                  // Set P(1, 0)
    sw         matrix,  +q0off1, v99                                  // Set Q(1, 0)

.Lphase1:
    move       cost,    COST_SUB                                      // Initialization value
    move       i,       1                                             // Initialization value
    add        mpp,     matrix,  +4                                   // Matrix line (index pp,i)
    add        mlp,     matrix,  +0                                   // Matrix line (index lp,i)
.Lphase1_iloop:
    sw         mpp,     +d0off,  cost                                 // Set D(pp, 0)
    move       mnScore, v99                                           // Initialization value

.Lphase1_v1:
    add        t0,      i,       -1                                   // Probably requires optimization here
    ror        t0,      t0,      2                                    //
    add        t1,      s1,      t0                                   //
    lbu        t1,      t1,      +0                                   //
    lsr        t0,      t0,      29                                   //
    lsr        v1,      t1,      t0                                   //

    move       j,       0                                             // Initialization value
    add        i,       i,       +(NB_DIAGS/2)                        // j loop bound
    lw         t0,      mpp,     +d0off                               // Get first D(pp, j-1)
.Lphase1_jloop:
    add        t0,      t0,      +COST_GAPO                           //
    lw         t1,      mpp,     +p0off                               // Get P(pp, j-1)
    add        t1,      t1,      +COST_GAPE                           //
    jltu       t0,      t1,      .Lphase1_min0                        // Min is D(pp, j-1)+COST_GAPO
    move       t0,      t1                                            // Min is P(pp, j-1)+COST_GAPE
.Lphase1_min0:
    sw         mpp,     +p1off,  t0                                   // Set P(pp, j)
    lw         QP,      mlp,     +d1off                               // Get D(lp, j)
    add        QP,      QP,      +COST_GAPO                           //
    lw         t1,      mlp,     +q1off                               // Get Q(lp, j)
    add        t1,      t1,      +COST_GAPE                           //
    jltu       QP,      t1,      .Lphase1_min1                        // Min is D(lp, j)+COST_GAPO
    move       QP,      t1                                            // Min is Q(lp, j)+COST_GAPE
.Lphase1_min1:
    sw         mpp,     +q1off,  QP                                   // Set Q(pp, j)
    jltu       QP,      t0,      .Lphase1_min2                        // Min is Q(pp, j)
    move       QP,      t0                                            // Min is P(pp, j)
.Lphase1_min2:

.Lphase1_v2:
    ror        t0,      j,       2                                    //
    add        t1,      s2,      t0                                   //
    lbu        t1,      t1,      +0                                   //
    lsr        t0,      t0,      29                                   //
    lsr        t1,      t1,      t0                                   //
    xor        t1,      t1,      v1                                   // base difference

.Lphase1_dcost:
    lw         t0,      mlp,     +d0off                               // Get D(lp, j-1)
    and        t1,      t1,      3, ?z, .Lphase1_nodiff               //
    add        t0,      t0,      +COST_SUB                            //
.Lphase1_nodiff:
    jltu       t0,      QP,      .Lphase1_min3                        //
    move       t0,      QP                                            //
.Lphase1_min3:
    sw         mpp,     +d1off,  t0                                   // Set D(pp, j), upcoming D(pp, j-1)
    jltu       mnScore, t0,      .Lphase1_min4                        // Min score saved
    move       mnScore, t0                                            // Min score updated
.Lphase1_min4:

.Lphase1_jnext:
    add        mpp,     mpp,     +LINE_SIZE                           // Loop increment
    add        mlp,     mlp,     +LINE_SIZE                           // Loop increment
    add        j,       j,       +1                                   // Loop increment
    jltu       j,       i,       .Lphase1_jloop                       // Next iteration

.Lphase1_earlyexit:
    sw         mpp,     +q1off,  v99                                  // Set Q(pp, j)
    sw         mpp,     +d1off,  v99                                  // Set D(pp, j)
    jltu       mxScore, mnScore, .Lrestore_registers                  // Early exit

.Lphase1_inext:
    add        cost,    cost,    +COST_SUB                            // Loop increment
    add        i,       i,       (1-NB_DIAGS/2)                       // Loop increment
    and        mpp,     i,       1, ?z, .Lphase1_odd                  //
    add        mpp,     matrix,  +4                                   // Switch vector
    add        mlp,     matrix,  +0, ?true, .Lphase1_itest            // Switch vector
.Lphase1_odd:
    add        mpp,     matrix,  +0                                   // Switch vector
    add        mlp,     matrix,  +4                                   // Switch vector
.Lphase1_itest:
    sub        t0,      i,       (NB_DIAGS/2+1), ?ltu, .Lphase1_iloop // Next iteration

.Lphase2:
.Lphase2_iloop:
    move       mnScore, v99                                           // Initialization value
    sw         mpp,     +p0off,  v99                                  // Set P(pp, j)
    sw         mpp,     +d0off,  v99                                  // Set D(pp, j)

.Lphase2_v1:
    add        t0,      i,       -1                                   // Probably requires optimization here
    ror        t0,      t0,      2                                    //
    add        t1,      s1,      t0                                   //
    lbu        t1,      t1,      +0                                   //
    lsr        t0,      t0,      29                                   //
    lsr        v1,      t1,      t0                                   //

    add        j,       i,       (-NB_DIAGS/2-1)                      // Initialization value
    add        i,       i,       (NB_DIAGS/2)                         // j loop bound
    lw         t0,      mpp,     +d0off                               // Get first D(pp, j-1)
.Lphase2_jloop:
    add        t0,      t0,      +COST_GAPO                           //
    lw         t1,      mpp,     +p0off                               // Get P(pp, j-1)
    add        t1,      t1,      +COST_GAPE                           //
    jltu       t0,      t1,      .Lphase2_min0                        // Min is D(pp, j-1)+COST_GAPO
    move       t0,      t1                                            // Min is P(pp, j-1)+COST_GAPE
.Lphase2_min0:
    sw         mpp,     +p1off,  t0                                   // Set P(pp, j)
    lw         QP,      mlp,     +d1off                               // Get D(lp, j)
    add        QP,      QP,      +COST_GAPO                           //
    lw         t1,      mlp,     +q1off                               // Get Q(lp, j)
    add        t1,      t1,      +COST_GAPE                           //
    jltu       QP,      t1,      .Lphase2_min1                        // Min is D(lp, j)+COST_GAPO
    move       QP,      t1                                            // Min is Q(lp, j)+COST_GAPE
.Lphase2_min1:
    sw         mpp,     +q1off,  QP                                   // Set Q(pp, j)
    jltu       QP,      t0,      .Lphase2_min2                        // Min is Q(pp, j)
    move       QP,      t0                                            // Min is P(pp, j)
.Lphase2_min2:

.Lphase2_v2:
    ror        t0,      j,       2                                    //
    add        t1,      s2,      t0                                   //
    lbu        t1,      t1,      +0                                   //
    lsr        t0,      t0,      29                                   //
    lsr        t1,      t1,      t0                                   //
    xor        t1,      t1,      v1                                   // base difference

.Lphase2_dcost:
    lw         t0,      mlp,     +d0off                               // Get D(lp, j-1)
    and        t1,      t1,      3, ?z, .Lphase2_nodiff               //
    add        t0,      t0,      +COST_SUB                            //
.Lphase2_nodiff:
    jltu       t0,      QP,      .Lphase2_min3                        //
    move       t0,      QP                                            //
.Lphase2_min3:
    sw         mpp,     +d1off,  t0                                   // Set D(pp, j), upcoming D(pp, j-1)
    jltu       mnScore, t0,      .Lphase2_min4                        // Min score saved
    move       mnScore, t0                                            // Min score updated
.Lphase2_min4:

.Lphase2_jnext:
    add        mpp,     mpp,     +LINE_SIZE                           // Loop increment
    add        mlp,     mlp,     +LINE_SIZE                           // Loop increment
    add        j,       j,       +1                                   // Loop increment
    jltu       j,       i,       .Lphase2_jloop                       // Next iteration

.Lphase2_earlyexit:
    sw         mpp,     +q1off,  v99                                  // Set Q(pp, j)
    sw         mpp,     +d1off,  v99                                  // Set D(pp, j)
    jltu       mxScore, mnScore, .Lrestore_registers                  // Early exit

.Lphase2_inext:
    add        mpp,     mpp,     (-LINE_SIZE*2*(NB_DIAGS/2))          // Loop increment
    add        mlp,     mlp,     (-LINE_SIZE*2*(NB_DIAGS/2))          // Loop increment
    xor        mpp,     mpp,     4                                    // Switch vector
    xor        mlp,     mlp,     4                                    // Switch vector
    add        i,       i,       (1-NB_DIAGS/2)                       // Loop increment
    add        t0,      len,     (-(NB_DIAGS/2))                      // Loop bound (to be optimized)
    jltu       i,       t0,      .Lphase2_iloop                       // Next iteration

.Lphase3:
    move       mnScore, v99                                           // Initialization value
    add        matrix,  mpp,     +LINE_SIZE                           // Next start

.Lphase3_iloop:
    sw         mpp,     +p0off,  v99                                  // Set P(pp, j)
    sw         mpp,     +d0off,  v99                                  // Set D(pp, j)

.Lphase3_v1:
    add        t0,      i,       -1                                   // Probably requires optimization here
    ror        t0,      t0,      2                                    //
    add        t1,      s1,      t0                                   //
    lbu        t1,      t1,      +0                                   //
    lsr        t0,      t0,      29                                   //
    lsr        v1,      t1,      t0                                   //

    add        j,       i,       (-NB_DIAGS/2-1)                      // Initialization value
    lw         t0,      mpp,     +d0off                               // Get first D(pp, j-1)
.Lphase3_jloop:
    add        t0,      t0,      +COST_GAPO                           //
    lw         t1,      mpp,     +p0off                               // Get P(pp, j-1)
    add        t1,      t1,      +COST_GAPE                           //
    jltu       t0,      t1,      .Lphase3_min0                        // Min is D(pp, j-1)+COST_GAPO
    move       t0,      t1                                            // Min is P(pp, j-1)+COST_GAPE
.Lphase3_min0:
    sw         mpp,     +p1off,  t0                                   // Set P(pp, j)
    lw         QP,      mlp,     +d1off                               // Get D(lp, j)
    add        QP,      QP,      +COST_GAPO                           //
    lw         t1,      mlp,     +q1off                               // Get Q(lp, j)
    add        t1,      t1,      +COST_GAPE                           //
    jltu       QP,      t1,      .Lphase3_min1                        // Min is D(lp, j)+COST_GAPO
    move       QP,      t1                                            // Min is Q(lp, j)+COST_GAPE
.Lphase3_min1:
    sw         mpp,     +q1off,  QP                                   // Set Q(pp, j)
    jltu       QP,      t0,      .Lphase3_min2                        // Min is Q(pp, j)
    move       QP,      t0                                            // Min is P(pp, j)
.Lphase3_min2:

.Lphase3_v2:
    ror        t0,      j,       2                                    //
    add        t1,      s2,      t0                                   //
    lbu        t1,      t1,      +0                                   //
    lsr        t0,      t0,      29                                   //
    lsr        t1,      t1,      t0                                   //
    xor        t1,      t1,      v1                                   // base difference

.Lphase3_dcost:
    lw         t0,      mlp,     +d0off                               // Get D(lp, j-1)
    and        t1,      t1,      3, ?z, .Lphase3_nodiff               //
    add        t0,      t0,      +COST_SUB                            //
.Lphase3_nodiff:
    jltu       t0,      QP,      .Lphase3_min3                        //
    move       t0,      QP                                            //
.Lphase3_min3:
    sw         mpp,     +d1off,  t0                                   // Set D(pp, j), upcoming D(pp, j-1)

.Lphase3_jnext:
    add        mpp,     mpp,     +LINE_SIZE                           // Loop increment
    add        mlp,     mlp,     +LINE_SIZE                           // Loop increment
    add        j,       j,       +1                                   // Loop increment
    jltu       j,       len,     .Lphase3_jloop                       // Next iteration

.Lphase3_min4:
    jltu       mnScore, t0,      .Lphase3_inext                       // Min is mnScore,
    move       mnScore, t0                                            // Min is D(pp, j-1)

.Lphase3_inext:
    sub        t0,      mpp,     matrix                               // Offset to next line
    sub        mlp,     mlp,     t0                                   // Loop increment
    xor        mpp,     matrix,  4                                    // Switch vector
    xor        mlp,     mlp,     4                                    // Switch vector
    add        matrix,  mpp,     +LINE_SIZE                           // Next start
    add        i,       i,       +1                                   // Loop increment
    jleu       i,       len,     .Lphase3_iloop                       // Next iteration

.Lphase4:
    xor        mpp,     mpp,     4                                    // Switch vector
    add        j,       len,     (-(NB_DIAGS/2))                      // Initialization value
.Lphase4_jloop:
    lw         t0,      mpp,     +d0off                               // Get D(pp,j)
    jltu       mnScore, t0,      .Lphase4_jnext                       // Min is mnScore
    move       mnScore, t0                                            // Min is D(pp,j)
.Lphase4_jnext:
    add        mpp,     mpp,     +LINE_SIZE                           // Loop increment
    add        j,       j,       +1                                   // Loop increment
    jleu       j,       len,     .Lphase4_jloop                       // Next iteration

.Lrestore_registers:
    ld         d14,     r22,     0
    jump       r23

.Lodpd_end:
.size odpd, .Lodpd_end-odpd
.section .stack_sizes,"o",@progbits,.text.odpd,unique,0
.long odpd
.byte 8
.section    .text.odpd,"ax",@progbits
