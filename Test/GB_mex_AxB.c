//------------------------------------------------------------------------------
// GB_mex_AxB: compute C=A*B, A'*B, A*B', or A'*B'
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2018, All Rights Reserved.
// http://suitesparse.com   See GraphBLAS/Doc/License.txt for license.

//------------------------------------------------------------------------------

// This is for testing only.  See GrB_mxm instead.  Returns a plain MATLAB
// matrix, in double.

#include "GB_mex.h"

#define USAGE "C = GB_mex_AxB (A, B, atranspose, btranspose, axb_method)"

#define FREE_ALL                        \
{                                       \
    GB_MATRIX_FREE (&A) ;               \
    GB_MATRIX_FREE (&Aconj) ;           \
    GB_MATRIX_FREE (&B) ;               \
    GB_MATRIX_FREE (&Bconj) ;           \
    GB_MATRIX_FREE (&C) ;               \
    GB_MATRIX_FREE (&Mask) ;            \
    GrB_free (&add) ;                   \
    GrB_free (&semiring) ;              \
    GB_mx_put_global (true) ;           \
}

//------------------------------------------------------------------------------

GrB_Info info ;
bool malloc_debug = false ;
bool ignore = false ;
bool atranspose = false ;
bool btranspose = false ;
GrB_Matrix A = NULL, B = NULL, C = NULL, Aconj = NULL, Bconj = NULL,
    Mask = NULL ;
GrB_Monoid add = NULL ;
GrB_Semiring semiring = NULL ;
int64_t anrows = 0 ;
int64_t ancols = 0 ;
int64_t bnrows = 0 ;
int64_t bncols = 0 ;

GrB_Desc_Value AxB_method = GxB_DEFAULT ;

//------------------------------------------------------------------------------

GrB_Info axb ( )
{

    // create the Semiring for regular z += x*y
    info = GrB_Monoid_new (&add, GrB_PLUS_FP64, (double) 0) ;
    if (info != GrB_SUCCESS) return (info) ;

    info = GrB_Semiring_new (&semiring, add, GrB_TIMES_FP64) ;
    if (info != GrB_SUCCESS)
    {
        GrB_free (&add) ;
        return (info) ;
    }

    /*
    double tic [2] ;
    simple_tic (tic) ;
    */

    // C = A*B, A'*B, A*B', or A'*B'
    info = GB_AxB_meta (&C, true /* CSC */,
        NULL /* no MT returned */,
        NULL /* no Mask */,
        A, B, semiring, /* GrB_PLUS_TIMES_FP64 */
        atranspose, btranspose, false, &ignore, AxB_method) ;

    /*
    double t = simple_toc (tic) ;
    printf ("GB_AxB_meta atrans %d btrans %d time %g method %d\n", atranspose, btranspose, t,
        GB_thread_local.AxB_method) ;
    */

    GrB_free (&add) ;
    GrB_free (&semiring) ;

    return (info) ;
}

//------------------------------------------------------------------------------

GrB_Info axb_complex ( )
{

    // C = A*B, complex case

    Aconj = NULL ;
    Bconj = NULL ;

    if (atranspose)
    {
        // Aconj = A
        info = GrB_Matrix_new (&Aconj, Complex, A->vlen, A->vdim) ;
        if (info != GrB_SUCCESS) return (info) ;
        info = GrB_apply (Aconj, NULL, NULL, Complex_conj, A, NULL) ;
        if (info != GrB_SUCCESS)
        {
            GrB_free (&Aconj) ;
            return (info) ;
        }
    }

    if (btranspose)
    {
        // Bconj = B
        info = GrB_Matrix_new (&Bconj, Complex, B->vlen, B->vdim) ;
        if (info != GrB_SUCCESS)
        {
            GrB_free (&Aconj) ;
            return (info) ;
        }

        info = GrB_apply (Bconj, NULL, NULL, Complex_conj, B, NULL) ;
        if (info != GrB_SUCCESS)
        {
            GrB_free (&Bconj) ;
            GrB_free (&Aconj) ;
            return (info) ;
        }

    }

    info = GB_AxB_meta (&C, true /*CSC*/,
        NULL /* no MT returned */,
        NULL /* no Mask */,
        (atranspose) ? Aconj : A,
        (btranspose) ? Bconj : B, Complex_plus_times,
        atranspose, btranspose, false, &ignore, AxB_method) ;

    GrB_free (&Bconj) ;
    GrB_free (&Aconj) ;

    return (info) ;
}

//------------------------------------------------------------------------------

void mexFunction
(
    int nargout,
    mxArray *pargout [ ],
    int nargin,
    const mxArray *pargin [ ]
)
{

    info = GrB_SUCCESS ;
    malloc_debug = GB_mx_get_global (true) ;
    ignore = false ;
    A = NULL ;
    B = NULL ;
    C = NULL ;
    Aconj = NULL ;
    Bconj = NULL ;
    Mask = NULL ;
    add = NULL ;
    semiring = NULL ;

    WHERE (USAGE) ;

    // check inputs
    if (nargout > 1 || nargin < 2 || nargin > 5)
    {
        mexErrMsgTxt ("Usage: " USAGE) ;
    }

    #define GET_DEEP_COPY ;
    #define FREE_DEEP_COPY ;

    if (mxIsComplex (pargin [0]))
    {
        // just for testing
        // METHOD (Complex_finalize ()) ;
        // METHOD (Complex_init ()) ;
    }

    // get A and B
    A = GB_mx_mxArray_to_Matrix (pargin [0], "A", false, true) ;
    B = GB_mx_mxArray_to_Matrix (pargin [1], "B", false, true) ;
    if (A == NULL || B == NULL)
    {
        FREE_ALL ;
        mexErrMsgTxt ("failed") ;
    }

    if (!A->is_csc || !B->is_csc)
    {
        mexErrMsgTxt ("A and B must be in CSC format") ;
    }

    // get the atranspose option
    GET_SCALAR (2, bool, atranspose, false) ;

    // get the btranspose option
    GET_SCALAR (3, bool, btranspose, false) ;

    // get the axb_method
    // 0 or not present: default
    // 1001: Gustavson
    // 1002: heap
    // 1003: dot
    GET_SCALAR (4, GrB_Desc_Value, AxB_method, GxB_DEFAULT) ;

    if (! ((AxB_method == GxB_DEFAULT) ||
        (AxB_method == GxB_AxB_GUSTAVSON) ||
        (AxB_method == GxB_AxB_HEAP) ||
        (AxB_method == GxB_AxB_DOT)))
    {
        mexErrMsgTxt ("unknown method") ;
    }

    // determine the dimensions
    anrows = (atranspose) ? NCOLS (A) : NROWS (A) ;
    ancols = (atranspose) ? NROWS (A) : NCOLS (A) ;
    bnrows = (btranspose) ? NCOLS (B) : NROWS (B) ;
    bncols = (btranspose) ? NROWS (B) : NCOLS (B) ;
    if (ancols != bnrows)
    {
        FREE_ALL ;
        mexErrMsgTxt ("invalid dimensions") ;
    }

    if (A->type == Complex)
    {
        METHOD (axb_complex ( )) ;
    }
    else
    {
        METHOD (axb ( )) ;
    }

    // return C to MATLAB
    pargout [0] = GB_mx_Matrix_to_mxArray (&C, "C AxB result", false) ;

    FREE_ALL ;
    GrB_finalize ( ) ;
}

