//------------------------------------------------------------------------------
// GB_mex_AdotB: compute C=spones(Mask).*(A'*B)
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2018, All Rights Reserved.
// http://suitesparse.com   See GraphBLAS/Doc/License.txt for license.

//------------------------------------------------------------------------------

// Returns a plain MATLAB sparse matrix, not a struct.  Only works in double.

#include "GB_mex.h"

#define USAGE "C = GB_mex_AdotB (A,B,Mask)"

#define FREE_ALL                        \
{                                       \
    GB_MATRIX_FREE (&A) ;               \
    GB_MATRIX_FREE (&Aconj) ;           \
    GB_MATRIX_FREE (&B) ;               \
    GB_MATRIX_FREE (&C) ;               \
    GB_MATRIX_FREE (&Mask) ;            \
    GrB_free (&add) ;                   \
    GrB_free (&semiring) ;              \
    GB_mx_put_global (true) ;           \
}

GrB_Matrix A = NULL, B = NULL, C = NULL, Aconj = NULL, Mask = NULL ;
GrB_Monoid add = NULL ;
GrB_Semiring semiring = NULL ;

//------------------------------------------------------------------------------

GrB_Info adotb_complex ( )
{
    GrB_Info info = GrB_Matrix_new (&Aconj, Complex, A->vlen, A->vdim) ;
    if (info != GrB_SUCCESS) return (info) ;
    info = GrB_apply (Aconj, NULL, NULL, Complex_conj, A, NULL) ;
    if (info != GrB_SUCCESS)
    {
        GrB_free (&Aconj) ;
        return (info) ;
    }
    info = GB_AxB_dot (&C, Mask, Aconj, B, Complex_plus_times,false) ;
    GrB_free (&Aconj) ;
    return (info) ;
}

//------------------------------------------------------------------------------

GrB_Info adotb ( ) 
{
    // create the Semiring for regular z += x*y
    GrB_Info info = GrB_Monoid_new (&add, GrB_PLUS_FP64, (double) 0) ;
    if (info != GrB_SUCCESS) return (info) ;
    info = GrB_Semiring_new (&semiring, add, GrB_TIMES_FP64) ;
    if (info != GrB_SUCCESS)
    {
        GrB_free (&add) ;
        return (info) ;
    }
    // C = A'*B
    info = GB_AxB_dot (&C, Mask, A, B,
        semiring /* GrB_PLUS_TIMES_FP64 */, false) ;
    GrB_free (&add) ;
    GrB_free (&semiring) ;
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

    bool malloc_debug = GB_mx_get_global (true) ;

    WHERE (USAGE) ;

    // check inputs
    if (nargout > 1 || nargin < 2 || nargin > 3)
    {
        mexErrMsgTxt ("Usage: " USAGE) ;
    }

    #define GET_DEEP_COPY ;
    #define FREE_DEEP_COPY ;

    GET_DEEP_COPY ;
    // get A and B (shallow copies)
    A = GB_mx_mxArray_to_Matrix (pargin [0], "A input", false, true) ;
    B = GB_mx_mxArray_to_Matrix (pargin [1], "B input", false, true) ;
    if (A == NULL)
    {
        FREE_ALL ;
        mexErrMsgTxt ("A failed") ;
    }
    if (B == NULL)
    {
        FREE_ALL ;
        mexErrMsgTxt ("B failed") ;
    }

    // get Mask (shallow copy)
    if (nargin > 2)
    {
        Mask = GB_mx_mxArray_to_Matrix (pargin [2], "Mask input", false, false);
    }

    if (A->vlen != B->vlen)
    {
        FREE_ALL ;
        mexErrMsgTxt ("inner dimensions of A'*B do not match") ;
    }

    if (A->type == Complex)
    {
        // C = A'*B, complex case
        METHOD (adotb_complex ( )) ;
        /*
        METHOD (GrB_Matrix_new (&Aconj, Complex, A->vlen, A->vdim)) ;
        METHOD (GrB_apply (Aconj, NULL, NULL, Complex_conj, A, NULL)) ;
        METHOD (GB_AxB_dot (&C, Mask, Aconj, B, Complex_plus_times,false));
        */
    }
    else
    {
        METHOD (adotb ( )) ;

        #if 0
        // create the Semiring for regular z += x*y
        METHOD (GrB_Monoid_new (&add, GrB_PLUS_FP64, (double) 0)) ;
        METHOD (GrB_Semiring_new (&semiring, add, GrB_TIMES_FP64)) ;
        // C = A'*B
        METHOD (GB_AxB_dot (&C, Mask, A, B, semiring
            /* GrB_PLUS_TIMES_FP64 */, false)) ;
        #endif

    }

    // return C to MATLAB
    pargout [0] = GB_mx_Matrix_to_mxArray (&C, "C AdotB result", false) ;

    FREE_ALL ;
}

