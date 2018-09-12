//------------------------------------------------------------------------------
// GB_kron: C<M> = accum (C, kron(A,B))
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2018, All Rights Reserved.
// http://suitesparse.com   See GraphBLAS/Doc/License.txt for license.

//------------------------------------------------------------------------------

// C<M> = accum (C, kron(A,B))

// The input matrices A and B are optionally transposed.

// Not user-callable.  Does the work for GxB_kron

#include "GB.h"

GrB_Info GB_kron                    // C<M> = accum (C, kron(A,B))
(
    GrB_Matrix C,                   // input/output matrix for results
    const bool C_replace,           // if true, clear C before writing to it
    const GrB_Matrix M,             // optional mask for C, unused if NULL
    const bool Mask_comp,           // if true, use ~M
    const GrB_BinaryOp accum,       // optional accum for Z=accum(C,T)
    const GrB_BinaryOp op,          // defines '*' for kron(A,B)
    const GrB_Matrix A,             // input matrix
    bool A_transpose,               // if true, use A' instead of A
    const GrB_Matrix B,             // input matrix
    bool B_transpose                // if true, use B' instead of B
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    ASSERT (ALIAS_OK3 (C, M, A, B)) ;

    RETURN_IF_NULL_OR_FAULTY (C) ;
    RETURN_IF_NULL_OR_FAULTY (A) ;
    RETURN_IF_NULL_OR_FAULTY (B) ;
    RETURN_IF_FAULTY (M) ;
    RETURN_IF_NULL_OR_FAULTY (op) ;
    RETURN_IF_FAULTY (accum) ;

    ASSERT_OK (GB_check (C, "C input for GB_kron", D0)) ;
    ASSERT_OK_OR_NULL (GB_check (M, "M for GB_kron", D0)) ;
    ASSERT_OK_OR_NULL (GB_check (accum, "accum for GB_kron", D0)) ;
    ASSERT_OK (GB_check (op, "op for GB_kron", D0)) ;
    ASSERT_OK (GB_check (A, "A for GB_kron", D0)) ;
    ASSERT_OK (GB_check (B, "B for GB_kron", D0)) ;

    // check domains and dimensions for C<M> = accum (C,T)
    GrB_Info info = GB_compatible (C->type, C, M, accum, op->ztype) ;
    if (info != GrB_SUCCESS)
    { 
        return (info) ;
    }

    // T=op(A,B) via op operator, so A and B must be compatible with z=op(a,b)
    info = GB_BinaryOp_compatible (op, NULL, A->type, B->type, 0) ;
    if (info != GrB_SUCCESS)
    { 
        return (info) ;
    }

    // delete any lingering zombies and assemble any pending tuples in A and B,
    // so that cnz = NNZ(A) * NNZ(B) can be computed.  Pending updates of C
    // and M are done after this check.
    WAIT (A) ;
    WAIT (B) ;

    // check the dimensions of C
    int64_t anrows = (A_transpose) ? NCOLS (A) : NROWS (A) ;
    int64_t ancols = (A_transpose) ? NROWS (A) : NCOLS (A) ;
    int64_t bnrows = (B_transpose) ? NCOLS (B) : NROWS (B) ;
    int64_t bncols = (B_transpose) ? NROWS (B) : NCOLS (B) ;
    GrB_Index cnrows, cncols, cnz = 0 ;
    bool ok = GB_Index_multiply (&cnrows, anrows,  bnrows) ;
    ok = ok && GB_Index_multiply (&cncols, ancols,  bncols) ;
    ok = ok && GB_Index_multiply (&cnz, NNZ (A), NNZ (B)) ;
    if (!ok || NROWS (C) != cnrows || NCOLS (C) != cncols)
    { 
        return (ERROR (GrB_DIMENSION_MISMATCH, (LOG, "%s:\n"
            "output is "GBd"-by-"GBd"; must be "GBd"-by-"GBd"\n"
            "first input is "GBd"-by-"GBd"%s with "GBd" entries\n"
            "second input is "GBd"-by-"GBd"%s with "GBd" entries",
            ok ? "Dimensions not compatible:" : "Problem too large:",
            NROWS (C), NCOLS (C), cnrows, cncols,
            anrows, ancols, A_transpose ? " (transposed)" : "", NNZ (A),
            bnrows, bncols, B_transpose ? " (transposed)" : "", NNZ (B)))) ;
    }

    // quick return if an empty mask is complemented
    RETURN_IF_QUICK_MASK (C, C_replace, M, Mask_comp) ;

    // delete any lingering zombies and assemble any pending tuples
    WAIT (C) ;
    WAIT (M) ;

    //--------------------------------------------------------------------------
    // transpose A and B if requested
    //--------------------------------------------------------------------------

    bool is_csc = C->is_csc ;
    if (is_csc != A->is_csc)
    { 
        // Flip the sense of A_transpose
        A_transpose = !A_transpose ;
    }
    if (is_csc != B->is_csc)
    { 
        // Flip the sense of B_transpose
        B_transpose = !B_transpose ;
    }

    GrB_Matrix AT = NULL ;
    if (A_transpose)
    {
        // AT = A' and typecast to op->xtype
        // transpose: shallow output ok, typecast, no op
        info = GB_transpose (&AT, op->xtype, is_csc, A, NULL) ;
        if (info != GrB_SUCCESS)
        { 
            return (info) ;
        }
        ASSERT_OK (GB_check (A , "A after AT kron", D0)) ;
        ASSERT_OK (GB_check (AT, "AT kron", D0)) ;
    }

    GrB_Matrix BT = NULL ;
    if (B_transpose)
    {
        // BT = B' and typecast to op->ytype
        // transpose: shallow output ok, typecast, no op
        info = GB_transpose (&BT, op->ytype, is_csc, B, NULL) ;
        if (info != GrB_SUCCESS)
        { 
            GB_MATRIX_FREE (&AT) ;
            return (info) ;
        }
        ASSERT_OK (GB_check (BT, "BT kron", D0)) ;
    }

    //--------------------------------------------------------------------------
    // T = kron(A,B)
    //--------------------------------------------------------------------------

    GrB_Matrix T ;
    info = GB_kron_kernel (&T, C->is_csc, op,
        A_transpose ? AT : A, B_transpose ? BT : B) ;

    // free workspace
    GB_MATRIX_FREE (&AT) ;
    GB_MATRIX_FREE (&BT) ;

    if (info != GrB_SUCCESS)
    { 
        return (info) ;
    }

    ASSERT_OK (GB_check (T, "T = kron(A,B)", D0)) ;

    //--------------------------------------------------------------------------
    // C<M> = accum (C,T): accumulate the results into C via the mask
    //--------------------------------------------------------------------------

    return (GB_accum_mask (C, M, NULL, accum, &T, C_replace, Mask_comp)) ;
}

