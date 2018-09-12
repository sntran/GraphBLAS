//------------------------------------------------------------------------------
// GrB_transpose: transpose a sparse matrix
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2018, All Rights Reserved.
// http://suitesparse.com   See GraphBLAS/Doc/License.txt for license.

//------------------------------------------------------------------------------

// C<M> = accum (C,A') or accum (C,A)

#include "GB.h"

GrB_Info GrB_transpose              // C<M> = accum(C,A') or accum(C,A)
(
    GrB_Matrix C,                   // input/output matrix for results
    const GrB_Matrix M,             // optional mask for C, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(C,T)
    const GrB_Matrix A,             // first input:  matrix A
    const GrB_Descriptor desc       // descriptor for C, M, and A
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    ASSERT (ALIAS_OK2 (C, M, A)) ;

    WHERE ("GrB_transpose (C, M, accum, A, desc)") ;
    RETURN_IF_NULL_OR_FAULTY (C) ;
    RETURN_IF_FAULTY (M) ;
    RETURN_IF_FAULTY (accum) ;
    RETURN_IF_NULL_OR_FAULTY (A) ;

    ASSERT_OK (GB_check (C, "C input for GrB_transpose", D0)) ;
    ASSERT_OK_OR_NULL (GB_check (M, "M for GrB_transpose", D0)) ;
    ASSERT_OK_OR_NULL (GB_check (accum, "accum for GrB_transpose", D0)) ;
    ASSERT_OK (GB_check (A, "A input for GrB_transpose", D0)) ;

    // get the descriptor
    GET_DESCRIPTOR (info, desc, C_replace, Mask_comp, A_transpose, xx1, xx2) ;

    // check domains and dimensions for C<M> = accum (C,T)
    info = GB_compatible (C->type, C, M, accum, A->type) ;
    if (info != GrB_SUCCESS)
    { 
        return (info) ;
    }

    // check the dimensions
    int64_t tnrows = (!A_transpose) ? NCOLS (A) : NROWS (A) ;
    int64_t tncols = (!A_transpose) ? NROWS (A) : NCOLS (A) ;
    if (NROWS (C) != tnrows || NCOLS (C) != tncols)
    { 
        return (ERROR (GrB_DIMENSION_MISMATCH, (LOG,
            "Dimensions not compatible:\n"
            "output is "GBd"-by-"GBd"\n"
            "input is "GBd"-by-"GBd"%s",
            NROWS (C), NCOLS (C),
            tnrows, tncols, (!A_transpose) ? " (transposed)" : ""))) ;
    }

    // quick return if an empty mask is complemented
    RETURN_IF_QUICK_MASK (C, C_replace, M, Mask_comp) ;

    // delete any lingering zombies and assemble any pending tuples
    WAIT (C) ;
    WAIT (M) ;
    WAIT (A) ;

    //--------------------------------------------------------------------------
    // T = A or A', where T can have the type of C or the type of A
    //--------------------------------------------------------------------------

    bool C_is_csc = C->is_csc ;
    if (C_is_csc != A->is_csc)
    { 
        // Flip the sense of A_transpose
        A_transpose = !A_transpose ;
    }

    GrB_Matrix T = NULL ;

    if (!A_transpose)
    {
        // T = A', the default behavior.  This step may seem counter-intuitive,
        // but method computes C<M>=A' by default when A_transpose is false.

        // Precasting:
        if (accum == NULL)
        { 
            // If there is no accum operator, T is transplanted into Z and
            // typecasted into the C->type during the transpose.
            // transpose: shallow output ok, typecast, no op
            info = GB_transpose (&T, C->type, C_is_csc, A, NULL) ;
        }
        else
        { 
            // If the accum operator is present, entries in the intersection of
            // T and C are typecasted into the accum->ytype, while entries in T
            // but not C are typecasted directly into C->type.  Thus, the
            // typecast of T (if any) must wait, and be done in call to
            // GB_add in GB_accum_mask.
            // transpose: shallow output ok, no typecast, no op
            info = GB_transpose (&T, A->type, C_is_csc, A, NULL) ;
        }

        // no operator; typecasting done if accum is NULL
    }
    else
    { 
        // T = A, a pure shallow copy; nothing at all is allocated.  No
        // typecasting is done since the types of T and A are the same.  If the
        // A_transpose descriptor is true, A is viewed as transposed first.
        // The method transposes A again, giving T=A''=A.  This is a double
        // transpose, so C<M>=A is computed, and no transpose is done.  T is
        // typecasted eventually, into the type of C if the types of T and C
        // differ.  That can be postponed at no cost since the following step
        // is free.
        info = GB_shallow_cast (&T, A->type, C_is_csc, A) ;
    }

    if (info != GrB_SUCCESS)
    { 
        ASSERT (T == NULL) ;
        return (info) ;
    }

    ASSERT (T->is_csc == C->is_csc) ;

    //--------------------------------------------------------------------------
    // C<M> = accum (C,T): accumulate the results into C via the mask M
    //--------------------------------------------------------------------------

    ASSERT_OK (GB_check (T, "T for GrB_transpose", D0)) ;
    return (GB_accum_mask (C, M, NULL, accum, &T, C_replace, Mask_comp)) ;
}

