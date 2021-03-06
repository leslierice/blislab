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
#include "bl_dgemm.h"

/*
 * --------------------------------------------------------------------------
 */

// C must be aligned
void bl_dgemm_str_naive(
        int    m,
        int    n,
        int    k,
        double *XAA,
        double *XAB,
        double *packA,
        int    lda,
        int    gammaA,
        double *XBA,
        double *XBB,
        double *packB,
        int    ldb,
        int    gammaB,
        double *CA,        // must be aligned
        double *CB,
        int    ldc,        // ldc must also be aligned
        int    gammaCA,
        int    gammaCB,
        int    bl_ic_nt
        )
{
    //malloc XA, XB, XC_tmp
    //printf( "m=%d,n=%d,k=%d,lda=%d,ldb=%d,ldc=%d\n", m, n, k, lda, ldb, ldc );

    double *XA    = (double*)malloc( sizeof(double) * lda * k );
    double *XB    = (double*)malloc( sizeof(double) * ldb * n );
    double *XC_tmp = bl_malloc_aligned( ldc, n + 4, sizeof(double) );

    memset( XC_tmp, 0, sizeof(double) * ldc * ( n + 4 ) );

    double dbl_gammaA, dbl_gammaB;
    double dbl_gammaCA, dbl_gammaCB;

    double *ptr_gammaA, *ptr_gammaB;
    double *ptr_gammaCA, *ptr_gammaCB;

    //XA = XAA + gammaA * XAB
    mkl_copym( m, k, XAA, XA, lda );
    if ( gammaA != 0 ) {
        dbl_gammaA = (double) gammaA;
        ptr_gammaA = &dbl_gammaA;
        mkl_axpym( m, k, ptr_gammaA, XAB, XA, lda);
    }

    //XB = XBA + gammaB * XBB
    mkl_copym( k, n, XBA, XB, ldb );
    if ( gammaB != 0 ) {
        dbl_gammaB = (double) gammaB;
        ptr_gammaB = &dbl_gammaB;
        mkl_axpym( k, n, ptr_gammaB, XBB, XB, ldb );
    }

    bl_dgemm( m, n, k, XA, lda, XB, ldb, XC_tmp, ldc );

    dbl_gammaCA = (double) gammaCA;
    ptr_gammaCA = &dbl_gammaCA;
    // XC_tmp: size: ldc * n, not m * n;
    mkl_axpym( m, n, ptr_gammaCA, XC_tmp, CA, ldc );

    if ( gammaCB != 0 ) {
        dbl_gammaCB = (double) gammaCB;
        ptr_gammaCB = &dbl_gammaCB;
        mkl_axpym( m, n, ptr_gammaCB, XC_tmp, CB, ldc );
    }

    free( XA     );
    free( XB     );
    free( XC_tmp );

}

/*
 * --------------------------------------------------------------------------
 */

void bl_dgemm_strassen_naive(
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

    // M1: c00 = 1*c00+1*(a00+a11)(b00+b11); c11 = 1*c11+1*(a00+a11)(b00+b11)
    bl_dgemm_str_naive(ms, ns, ks, &XA[ 0 ], &XA[ ms + lda * ks ], packA, lda, 1, &XB[ 0 ], &XB[ ks + ldb * ns ], packB, ldb, 1, &C[ 0 ], &C[ ms + ldc * ns ], ldc, 1, 1, bl_ic_nt);
    // M2: c10 = 1*c10+1*(a10+a11)b00; c11 = 1*c11-1*(a10+a11)b00
    bl_dgemm_str_naive(ms, ns, ks, &XA[ ms ], &XA[ ms + lda * ks ], packA, lda, 1, &XB[ 0 ], &XB[ 0 ], packB, ldb, 0, &C[ ms ], &C[ ms + ldc * ns ], ldc, 1, -1, bl_ic_nt);
    // M3: c01 = 1*c01+1*a00(b01-b11); c11 = 1*c11+1*a00(b01-b11)
    bl_dgemm_str_naive(ms, ns, ks, &XA[ 0 ], &XA[ 0 ], packA, lda, 0, &XB[ ldb * ns ], &XB[ ks + ldb * ns ], packB, ldb, -1, &C[ ldc * ns ], &C[ ms + ldc * ns ], ldc, 1, 1, bl_ic_nt);
    // M4: c00 = 1*c00+1*a11(b10-b00); c10 = 1*c10+1*a11(b10-b00)
    bl_dgemm_str_naive(ms, ns, ks, &XA[ ms + lda * ks ], &XA[ 0 ], packA, lda, 0, &XB[ ks ], &XB[ 0 ], packB, ldb, -1, &C[ 0 ], &C[ ms ], ldc, 1, 1, bl_ic_nt);
    // M5: c00 = 1*c00-1*(a00+a01)b11; c01 = 1*c01+1*(a00+a01)b11
    bl_dgemm_str_naive(ms, ns, ks, &XA[ 0 ], &XA[ lda * ks ], packA, lda, 1, &XB[ ks + ldb * ns ], &XB[ 0 ], packB, ldb, 0, &C[ ldc * ns ], &C[ 0 ], ldc, 1, -1, bl_ic_nt);
    // M6: c11 = 1*c11+(a10-a00)(b00+b01)
    bl_dgemm_str_naive(ms, ns, ks, &XA[ ms ], &XA[ 0 ], packA, lda, -1, &XB[ 0 ], &XB[ ldb * ns ], packB, ldb, 1, &C[ ms + ldc * ns ], &C[ 0 ], ldc, 1, 0, bl_ic_nt);
    // M7: c00 = 1*c00+(a01-a11)(b10+b11)
    bl_dgemm_str_naive(ms, ns, ks, &XA[ lda * ks ], &XA[ ms + lda * ks ], packA, lda, -1, &XB[ ks ], &XB[ ks + ldb * ns ], packB, ldb, 1, &C[ 0 ], &C[ 0 ], ldc, 1, 0, bl_ic_nt);

    free( packA );
    free( packB );
}


//bf = m_R = 8
