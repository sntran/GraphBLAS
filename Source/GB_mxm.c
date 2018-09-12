//------------------------------------------------------------------------------
// GB_mxm: matrix-matrix multiply for GrB_mxm, GrB_mxv, and GrB_vxm
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2018, All Rights Reserved.
// http://suitesparse.com   See GraphBLAS/Doc/License.txt for license.

//------------------------------------------------------------------------------

// C<M> = accum (C,A*B) and variations.

// This function is not user-callable.  It does the work for user-callable
// functions GrB_mxm, GrB_mxv, and GrB_vxm.

#include "GB.h"

GrB_Info GB_mxm                     // C<M> = A*B
(
    GrB_Matrix C,                   // input/output matrix for results
    const bool C_replace,           // if true, clear C before writing to it
    const GrB_Matrix M,             // optional mask for C, unused if NULL
    const bool Mask_comp,           // if true, use ~M
    const GrB_BinaryOp accum,       // optional accum for Z=accum(C,T)
    const GrB_Semiring semiring,    // defines '+' and '*' for C=A*B
    const GrB_Matrix A,             // input matrix
    const bool A_transpose,         // if true, use A' instead of A
    const GrB_Matrix B,             // input matrix
    const bool B_transpose,         // if true, use B' instead of B
    const bool flipxy,              // if true, do z=fmult(b,a) vs fmult(a,b)
    const GrB_Desc_Value AxB_method // for auto vs user selection of methods
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    ASSERT (ALIAS_OK3 (C, M, A, B)) ;

    RETURN_IF_FAULTY (accum) ;
    RETURN_IF_NULL_OR_FAULTY (semiring) ;

    ASSERT_OK (GB_check (C, "C input for GB_mxm", D0)) ;
    ASSERT_OK_OR_NULL (GB_check (M, "M for GB_mxm", D0)) ;
    ASSERT_OK_OR_NULL (GB_check (accum, "accum for GB_mxm", D0)) ;
    ASSERT_OK (GB_check (semiring, "semiring for GB_mxm", D0)) ;
    ASSERT_OK (GB_check (A, "A for GB_mxm", D0)) ;
    ASSERT_OK (GB_check (B, "B for GB_mxm", D0)) ;

    // check domains and dimensions for C<M> = accum (C,T)
    GrB_Type T_type = semiring->add->op->ztype ;
    GrB_Info info = GB_compatible (C->type, C, M, accum, T_type) ;
    if (info != GrB_SUCCESS)
    { 
        return (info) ;
    }

    // T=A*B via semiring, so A and B must be compatible with semiring->mult
    if (flipxy)
    { 
        // z=fmult(b,a), for entries a from A, and b from B
        info = GB_BinaryOp_compatible (semiring->multiply,
                                        NULL, B->type, A->type, 0) ;
    }
    else
    { 
        // z=fmult(a,b), for entries a from A, and b from B
        info = GB_BinaryOp_compatible (semiring->multiply,
                                        NULL, A->type, B->type, 0) ;
    }
    if (info != GrB_SUCCESS)
    { 
        return (info) ;
    }

    // check the dimensions
    int64_t anrows = (A_transpose) ? NCOLS (A) : NROWS (A) ;
    int64_t ancols = (A_transpose) ? NROWS (A) : NCOLS (A) ;
    int64_t bnrows = (B_transpose) ? NCOLS (B) : NROWS (B) ;
    int64_t bncols = (B_transpose) ? NROWS (B) : NCOLS (B) ;
    if (ancols != bnrows || NROWS (C) != anrows || NCOLS (C) != bncols)
    {
        return (ERROR (GrB_DIMENSION_MISMATCH, (LOG,
            "Dimensions not compatible:\n"
            "output is "GBd"-by-"GBd"\n"
            "first input is "GBd"-by-"GBd"%s\n"
            "second input is "GBd"-by-"GBd"%s",
            NROWS (C), NCOLS (C),
            anrows, ancols, A_transpose ? " (transposed)" : "",
            bnrows, bncols, B_transpose ? " (transposed)" : ""))) ;
    }

    // quick return if an empty mask is complemented
    RETURN_IF_QUICK_MASK (C, C_replace, M, Mask_comp) ;

    // delete any lingering zombies and assemble any pending tuples
    WAIT (C) ;
    WAIT (M) ;
    WAIT (A) ;
    WAIT (B) ;

    //--------------------------------------------------------------------------
    // T = A*B, A'*B, A*B', or A'*B', also using the mask to cut time and memory
    //--------------------------------------------------------------------------

    // the mask is used to cut time and memory usage for GB_AxB_meta,
    // but only if it is not complemented.
    bool mask_applied = false ;
    bool C_is_csc = C->is_csc ;
    GrB_Matrix T, MT = NULL ;
    info = GB_AxB_meta (&T, C_is_csc, &MT, ((Mask_comp) ? NULL : M), A, B,
        semiring, A_transpose, B_transpose, flipxy, &mask_applied, AxB_method) ;
    if (info != GrB_SUCCESS)
    { 
        // out of memory
        ASSERT (T == NULL) ;
        ASSERT (MT == NULL) ;
        return (info) ;
    }

    ASSERT_OK (GB_check (T, "T=A*B from GB_AxB_meta", D0)) ;
    ASSERT_OK_OR_NULL (GB_check (MT, "MT from GB_AxB_meta", D0)) ;
    ASSERT (!ZOMBIES (T)) ;
    ASSERT (!PENDING (T)) ;

    //--------------------------------------------------------------------------
    // C<M> = accum (C,T): accumulate the results into C via the mask
    //--------------------------------------------------------------------------

    if ((accum == NULL) && (C->is_csc == T->is_csc)
        && (M == NULL || (M != NULL && mask_applied))
        && (C_replace || NNZ (C) == 0))
    { 
        // C = 0 ; C = (ctype) T ; with the same CSR/CSC format.
        // The mask M (if any) has already been applied in GB_AxB_meta.
        // If C is also empty, or to be cleared anyway, and if accum is not
        // present, then T can be transplanted directly into C, as C = (ctype)
        // T, typecasting if needed.  If no typecasting is done then this takes
        // no time at all and is a pure transplant.  Also conform C to its
        // desired hypersparsity.
        GB_MATRIX_FREE (&MT) ;
        return (GB_transplant_conform (C, C->type, &T)) ;
    }
    else
    { 
        // C<M> = accum (C,T)
        // GB_accum_mask also conforms C to its desired hypersparsity
        info = GB_accum_mask (C, M, MT, accum, &T, C_replace, Mask_comp) ;
        GB_MATRIX_FREE (&MT) ;
        return (info) ;
    }
}

