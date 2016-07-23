/*
 * --------------------------------------------------------------------------
 * BLISLAB
 * --------------------------------------------------------------------------
 * Copyright (C) 2016, The University of Texas at Austin
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *  - Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *  - Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *  - Neither the name of The University of Texas nor the names of its
 *    contributors may be used to endorse or promote products derived
 *    from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *
 * bl_dgemm.c
 *
 *
 * Purpose:
 * this is the main file of blislab dgemm.
 *
 * Todo:
 *
 *
 * Modification:
 *
 *
 * */

#include <stdio.h>
#include <omp.h>

#include "bl_dgemm_kernel.h"
#include "bl_dsymm.h"

inline void packA_add_mcxkc_d_str_ab(
        int    m,
        int    k,
        double *XAA,
        int    strucAA,
        double *XAB,
        int    strucAB,
        double gamma,
        int    ldXA,
        int    inc_aa,
        int    inc_ab,
        int    ic_aa,
        int    ic_ab,
        int    pc_aa,
        int    pc_ab,
        int    offseta,
        double *packA
        )
{
    int    i, p;
    double *aa_pntr[ DGEMM_MR ];
    double *ab_pntr[ DGEMM_MR ];

    int diff_aa = ic_aa - pc_aa;
    int off_aa = diff_aa - offseta + 1;
    int diff_ab = ic_ab - pc_ab;
    int off_ab = diff_ab - offseta + 1;

    if ( strucAA == BL_TRANS_SYMM ) {
        for ( i = 0; i < min( m, off_aa ); i ++ ) {
            aa_pntr[ i ] = XAA + ldXA * ( offseta + i );
        }

        for ( i = min( m, off_aa ); i < m; i ++ ) {
            aa_pntr[ i ] = XAA + inc_aa * ( offseta + i - diff_aa ) + ldXA * diff_aa;
        }
    }
    else {
        for ( i = 0; i < m; i ++ ) {
            aa_pntr[ i ] = XAA + inc_aa * ( offseta + i );
        }

        for ( i = m; i < DGEMM_MR; i ++ ) {
            aa_pntr[ i ] = XAA + inc_aa * ( offseta + 0 );
        }
    }
    if ( strucAB == BL_TRANS_SYMM ) {
        for ( i = 0; i < min( m, off_ab ); i ++ ) {
            ab_pntr[ i ] = XAB + ldXA * ( offseta + i );
        }

        for ( i = min( m, off_ab ); i < m; i ++ ) {
            ab_pntr[ i ] = XAB + inc_ab * ( offseta + i - diff_ab ) + ldXA * diff_ab;
        }
    }
    else {
        for ( i = 0; i < m; i ++ ) {
            ab_pntr[ i ] = XAB + inc_ab * ( offseta + i );
        }

        for ( i = m; i < DGEMM_MR; i ++ ) {
            ab_pntr[ i ] = XAB + inc_ab * ( offseta + 0 );
        }
    }

    for ( p = 0; p < k; p ++ ) {
        for ( i = 0; i < DGEMM_MR; i ++ ) {
            *packA = *aa_pntr[ i ] + gamma * *ab_pntr[ i ];
            packA ++;

            if ( strucAA == BL_TRANS || ( strucAA == BL_SYMM && ( ic_aa + offseta + i <= pc_aa + p ) ) || ( strucAA == BL_TRANS_SYMM && i - p <= diff_aa ) ) {
                aa_pntr[ i ] ++;
            }
            else {
                aa_pntr[ i ] = aa_pntr[ i ] + ldXA;
            }
            if ( strucAB == BL_TRANS || ( strucAB == BL_SYMM && ( ic_ab + offseta + i <= pc_ab + p ) ) || ( strucAB == BL_TRANS_SYMM && i - p <= diff_ab ) ) {
                ab_pntr[ i ] ++;
            }
            else {
                ab_pntr[ i ] = ab_pntr[ i ] + ldXA;
            }
        }
    }
}


/*
 * --------------------------------------------------------------------------
 */

inline void packA_mcxkc_d_str_ab(
        int    m,
        int    k,
        double *XA,
        int    struc,
        int    ldXA,
        int    inc,
        int    ic,
        int    pc,
        int    offseta,
        double *packA
        )
 {
    int    i, p;
    double *a_pntr[ DGEMM_MR ];

    int diff = ic - pc;
    int off = diff - offseta + 1;

    if ( struc == BL_TRANS_SYMM ) {
        for ( i = 0; i < min( m, off ); i ++ ) {
            a_pntr[ i ] = XA + ldXA * ( offseta + i );
        }

        for ( i = min( m, off ); i < m; i ++ ) {
            a_pntr[ i ] = XA + inc * ( offseta + i - diff ) + ldXA * diff;
        }
    }
    else {
        for ( i = 0; i < m; i ++ ) {
            a_pntr[ i ] = XA + inc * ( offseta + i );
        }

        for ( i = m; i < DGEMM_MR; i ++ ) {
            a_pntr[ i ] = XA + inc * ( offseta + 0 );
        }
    }

    for ( p = 0; p < k; p ++ ) {
        for ( i = 0; i < DGEMM_MR; i ++ ) {
            *packA = *a_pntr[ i ];
            packA ++;

            if ( struc == BL_TRANS || ( struc == BL_SYMM && ( ic + offseta + i <= pc + p ) ) || ( struc == BL_TRANS_SYMM && i - p <= diff ) ) {
                a_pntr[ i ] = a_pntr[ i ] ++;
            }
            else {
                a_pntr[ i ] = a_pntr[ i ] + ldXA;
            }
        }
    }
}

/*
 * --------------------------------------------------------------------------
 */

inline void packB_add_kcxnc_d_str_ab(
        int    n,
        int    k,
        double *XBA,
        double *XBB,
        double gamma,
        int    ldXB, // ldXB is the original k
        int    offsetb,
        double *packB
        )
{
    int    j, p;
    double *ba_pntr[ DGEMM_NR ];
    double *bb_pntr[ DGEMM_NR ];

    for ( j = 0; j < n; j ++ ) {
        ba_pntr[ j ] = XBA + ldXB * ( offsetb + j );
        bb_pntr[ j ] = XBB + ldXB * ( offsetb + j );
    }

    for ( j = n; j < DGEMM_NR; j ++ ) {
        ba_pntr[ j ] = XBA + ldXB * ( offsetb + 0 );
        bb_pntr[ j ] = XBB + ldXB * ( offsetb + 0 );
    }

    for ( p = 0; p < k; p ++ ) {
        for ( j = 0; j < DGEMM_NR; j ++ ) {
            *packB ++ = (*ba_pntr[ j ] ++) + gamma * (*bb_pntr[ j ] ++);
        }
    }
}

/*
 * --------------------------------------------------------------------------
 */

inline void packB_kcxnc_d_str_ab(
        int    n,
        int    k,
        double *XB,
        int    ldXB, // ldXB is the original k
        int    offsetb,
        double *packB
        )
{
    int    j, p;
    double *b_pntr[ DGEMM_NR ];

    for ( j = 0; j < n; j ++ ) {
        b_pntr[ j ] = XB + ldXB * ( offsetb + j );
    }

    for ( j = n; j < DGEMM_NR; j ++ ) {
        b_pntr[ j ] = XB + ldXB * ( offsetb + 0 );
    }

    for ( p = 0; p < k; p ++ ) {
        for ( j = 0; j < DGEMM_NR; j ++ ) {
            *packB ++ = *b_pntr[ j ] ++;
        }
    }
}

/*
 * --------------------------------------------------------------------------
 */

void bl_macro_kernel_str_ab(
        int    m,
        int    n,
        int    k,
        double *packA,
        double *packB,
        double *CA,
        double *CB,
        int    ldc,
        double gammaCA,
        double gammaCB
        )
{
    int bl_ic_nt;
    int    i, ii, j;
    aux_t  aux;
    char *str;

    aux.b_next = packB;

    // We can also parallelize with OMP here.
    //// sequential is the default situation
    //bl_ic_nt = 1;
    //// check the environment variable
    //str = getenv( "BLISLAB_IC_NT" );
    //if ( str != NULL ) {
    //    bl_ic_nt = (int)strtol( str, NULL, 10 );
    //}
    //#pragma omp parallel for num_threads( bl_ic_nt ) private( j, i, aux )
    for ( j = 0; j < n; j += DGEMM_NR ) {                        // 2-th loop around micro-kernel
        aux.n  = min( n - j, DGEMM_NR );
        for ( i = 0; i < m; i += DGEMM_MR ) {                    // 1-th loop around micro-kernel
            aux.m = min( m - i, DGEMM_MR );
            if ( i + DGEMM_MR >= m ) {
                aux.b_next += DGEMM_NR * k;
            }

            double *c_tmp = calloc( DGEMM_MR * DGEMM_NR * k, sizeof(double) );

            ( *bl_micro_kernel_strassen ) (
                    k,
                    &packA[ i * k ],
                    &packB[ j * k ],
                    c_tmp,
                    CB,
                    (unsigned long long) DGEMM_MR,
                    gammaCA,
                    0,
                    &aux
                    );

            mkl_axpym( DGEMM_MR, DGEMM_NR, &gammaCA, c_tmp, DGEMM_MR, &CA[ j * ldc + i ], ldc );
            if ( gammaCB != 0 ) {
                mkl_axpym( DGEMM_MR, DGEMM_NR, &gammaCB, c_tmp, DGEMM_MR, &CB[ j * ldc + i ], ldc );
            }

            free( c_tmp );

        }                                                        // 1-th loop around micro-kernel
    }                                                            // 2-th loop around micro-kernel
}

/*
 * --------------------------------------------------------------------------
 */

// C must be aligned
void bl_dsymm_str_ab(
        int    m,
        int    n,
        int    k,
        double *XAA,
        int    strucAA,
        double *XAB,
        int    strucAB,
        double *packA,
        int    lda,
        double gammaA,
        double *XBA,
        double *XBB,
        double *packB,
        int    ldb,
        double gammaB,
        double *CA,        // must be aligned
        double *CB,
        int    ldc,        // ldc must also be aligned
        double gammaCA,
        double gammaCB,
        int    bl_ic_nt
        )
{
    int    i, j, p;
    int    ic, ib, jc, jb, pc, pb;
    int    ic_aa, ic_ab, pc_aa, pc_ab;
    int    ir, jr;
    int    inc_aa, inc_ab;
    int    strucAA_p, strucAB_p;

    for ( jc = 0; jc < n; jc += DGEMM_NC ) {                                       // 5-th loop around micro-kernel
        jb = min( n - jc, DGEMM_NC );
        for ( pc = 0; pc < k; pc += DGEMM_KC ) {                                   // 4-th loop around micro-kernel
            pb = min( k - pc, DGEMM_KC );

            #pragma omp parallel for num_threads( bl_ic_nt ) private( jr )
            for ( j = 0; j < jb; j += DGEMM_NR ) {
                if ( gammaB == 0 ) {
                    packB_kcxnc_d_str_ab(
                            min( jb - j, DGEMM_NR ),
                            pb,
                            &XBA[ pc ],
                            ldb,
                            jc + j,
                            &packB[ j * pb ]
                            );
                }
                else {
                    packB_add_kcxnc_d_str_ab(
                            min( jb - j, DGEMM_NR ),
                            pb,
                            &XBA[ pc ],
                            &XBB[ pc ],
                            gammaB,
                            ldb,
                            jc + j,
                            &packB[ j * pb ]
                            );
                }
            }

            //#pragma omp parallel for num_threads( bl_ic_nt ) private( ic, ib, i, ir )
            #pragma omp parallel num_threads( bl_ic_nt ) private( ic, ib, i, ir, ic_aa, pc_aa, inc_aa, strucAA_p, ic_ab, pc_ab, inc_ab, strucAB_p )
            {
                int     tid      = omp_get_thread_num();
                int     my_start;
                int     my_end;

                bl_get_range( m, DGEMM_MR, &my_start, &my_end );

                for ( ic = my_start; ic < my_end; ic += DGEMM_MC ) {              // 3-rd loop around micro-kernel

                    ib = min( my_end - ic, DGEMM_MC );

                    ic_aa = ic;
                    pc_aa = pc;
                    inc_aa = 1;
                    strucAA_p = strucAA;

                    if ( strucAA == BL_TRANS ) {
                        inc_aa = lda;
                        if ( pc != ic ) {
                            ic_aa = pc;
                            pc_aa = ic;
                        }
                    }
                    else if ( strucAA == BL_SYMM && pc > ic ) {
                        ic_aa = pc;
                        pc_aa = ic;
                        if (ic_aa - pc_aa <= DGEMM_MC) {
                            strucAA_p = BL_TRANS_SYMM;
                        }
                        else{
                            strucAA_p = BL_TRANS;
                            inc_aa = lda;
                        }

                    }

                    ic_ab = ic;
                    pc_ab = pc;
                    inc_ab = 1;
                    strucAB_p = strucAB;

                    if ( strucAB == BL_TRANS ) {
                        inc_ab = lda;
                        if ( pc != ic ) {
                            ic_ab = pc;
                            pc_ab = ic;
                        }
                    }
                    else if ( strucAB == BL_SYMM && pc > ic ) {
                        ic_ab = pc;
                        pc_ab = ic;
                        if (ic_ab - pc_ab <= DGEMM_MC) {
                            strucAB_p = BL_TRANS_SYMM;
                        }
                        else{
                            strucAB_p = BL_TRANS;
                            inc_ab = lda;
                        }
                    }

                    for ( i = 0; i < ib; i += DGEMM_MR ) {
                        if ( strucAA_p == BL_TRANS_SYMM && i >= ic_aa - pc_aa )  {
                            strucAA_p = BL_SYMM;
                            ic_aa = ic;
                            pc_aa = pc;
                        }
                        if ( strucAB_p == BL_TRANS_SYMM && i >= ic_ab - pc_ab )  {
                            strucAB_p = BL_SYMM;
                            ic_ab = ic;
                            pc_ab = pc;
                        }
                        if ( gammaA == 0 ) {
                            packA_mcxkc_d_str_ab(
                                    min( ib - i, DGEMM_MR ),
                                    pb,
                                    &XAA[ pc_aa * lda + ic_aa ],
                                    strucAA_p,
                                    lda,
                                    inc_aa,
                                    ic_aa,
                                    pc_aa,
                                    i,
                                    &packA[ tid * DGEMM_MC * pb + i * pb ]
                                    );
                        }
                        else {
                            packA_add_mcxkc_d_str_ab(
                                    min( ib - i, DGEMM_MR ),
                                    pb,
                                    &XAA[ pc_aa * lda + ic_aa ],
                                    strucAA_p,
                                    &XAB[ pc_ab * lda + ic_ab ],
                                    strucAB_p,
                                    gammaA,
                                    lda,
                                    inc_aa,
                                    inc_ab,
                                    ic_aa,
                                    ic_ab,
                                    pc_aa,
                                    pc_ab,
                                    i,
                                    &packA[ tid * DGEMM_MC * pb + i * pb ]
                                    );
                        }
                    }

                    bl_macro_kernel_str_ab(
                            ib,
                            jb,
                            pb,
                            packA  + tid * DGEMM_MC * pb,
                            packB,
                            &CA[ jc * ldc + ic ],
                            &CB[ jc * ldc + ic ],
                            ldc,
                            gammaCA,
                            gammaCB
                            );

                }                                                                // End 3.rd loop around micro-kernel

            }
        }                                                                        // End 4.th loop around micro-kernel
    }                                                                            // End 5.th loop around micro-kernel
}

/*
 * --------------------------------------------------------------------------
 */

void bl_dsymm_strassen_ab(
        int m,
        int n,
        int k,
        double *XA,
        int lda,
        double *XB,
        int ldb,
        double *C,
        int ldc
        )
{
    double *packA, *packB;
    char   *str;
    int bl_ic_nt;

    // Early return if possible
    if ( m == 0 || n == 0 || k == 0 ) {
        printf( "bl_dgemm(): early return\n" );
        return;
    }

    int	  ms=m/2, ns=n/2, ks=k/2;

    // sequential is the default situation
    bl_ic_nt = 1;
    // check the environment variable
    str = getenv( "BLISLAB_IC_NT" );
    if ( str != NULL ) {
        bl_ic_nt = (int)strtol( str, NULL, 10 );
    }

    // Allocate packing buffers
    packA  = bl_malloc_aligned( DGEMM_KC, ( DGEMM_MC + 1 ) * bl_ic_nt, sizeof(double) );
    packB  = bl_malloc_aligned( DGEMM_KC, ( DGEMM_NC + 1 )           , sizeof(double) );

    bl_dsymm_str_ab(ms, ns, ks, &XA[ 0 ], BL_SYMM, &XA[ ms + lda * ks ], BL_SYMM, packA, lda, 1, &XB[ 0 ], &XB[ ks + ldb * ns ], packB, ldb, 1, &C[ 0 ], &C[ ms + ldc * ns ], ldc, 1, 1, bl_ic_nt);
    bl_dsymm_str_ab(ms, ns, ks, &XA[ ms ], BL_GEN, &XA[ ms + lda * ks ], BL_SYMM, packA, lda, 1, &XB[ 0 ], &XB[ 0 ], packB, ldb, 0, &C[ ms ], &C[ ms + ldc * ns ], ldc, 1, -1, bl_ic_nt);
    bl_dsymm_str_ab(ms, ns, ks, &XA[ 0 ], BL_SYMM, &XA[ 0 ], BL_SYMM, packA, lda, 0, &XB[ ldb * ns ], &XB[ ks + ldb * ns ], packB, ldb, -1, &C[ ldc * ns ], &C[ ms + ldc * ns ], ldc, 1, 1, bl_ic_nt);
    bl_dsymm_str_ab(ms, ns, ks, &XA[ ms + lda * ns ], BL_SYMM, &XA[ 0 ], BL_SYMM, packA, lda, 0, &XB[ ks ], &XB[ 0 ], packB, ldb, -1, &C[ 0 ], &C[ ms ], ldc, 1, 1, bl_ic_nt);
    bl_dsymm_str_ab(ms, ns, ks, &XA[ 0 ], BL_SYMM, &XA[ ms ], BL_TRANS, packA, lda, 1, &XB[ ks + ldb * ns ], &XB[ 0 ], packB, ldb, 0, &C[ ldc * ns ], &C[ 0 ], ldc, 1, -1, bl_ic_nt);
    bl_dsymm_str_ab(ms, ns, ks, &XA[ ms ], BL_GEN, &XA[ 0 ], BL_SYMM, packA, lda, -1, &XB[ 0 ], &XB[ ldb * ns ], packB, ldb, 1, &C[ ms + ldc * ns ], &C[ 0 ], ldc, 1, 0, bl_ic_nt);
    bl_dsymm_str_ab(ms, ns, ks, &XA[ ms ], BL_TRANS, &XA[ ms + lda * ks ], BL_SYMM, packA, lda, -1, &XB[ ks ], &XB[ ks + ldb * ns ], packB, ldb, 1, &C[ 0 ], &C[ 0 ], ldc, 1, 0, bl_ic_nt);

    free( packA );
    free( packB );
}


//bf = m_R = 8
