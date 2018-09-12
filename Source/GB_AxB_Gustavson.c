//------------------------------------------------------------------------------
// GB_AxB_Gustavson: C=A*B or C<M>=A*B, gather/scatter-based saxpy method
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2018, All Rights Reserved.
// http://suitesparse.com   See GraphBLAS/Doc/License.txt for license.

//------------------------------------------------------------------------------

// This method is agnostic to the CSR/CSC format.  The format of C is set
// to CSC but this is a placeholder that will be changed in GB_AxB_meta.

// This function is called only by GB_AxB_meta.

#include "GB.h"
#ifndef GBCOMPACT
#include "GB_heap.h"
#include "GB_AxB__semirings.h"
#endif

// C=A*B is successful, just free temporary matrices
#define FREE_WORK                           \
{                                           \
    GB_MATRIX_FREE (&A2) ;                  \
    GB_MATRIX_FREE (&B2) ;                  \
}

// C=A*B failed, free everything, including C, and all thread-local workspace
#define FREE_ALL                            \
{                                           \
    FREE_WORK ;                             \
    GB_MATRIX_FREE (Chandle) ;              \
    GB_wfree ( ) ;                          \
}

#define OK(method)                          \
{                                           \
    info = method ;                         \
    if (info != GrB_SUCCESS)                \
    {                                       \
        FREE_ALL ;                          \
        return (info) ;                     \
    }                                       \
}

GrB_Info GB_AxB_Gustavson           // C=A*B or C<M>=A*B, Gustavson's method
(
    GrB_Matrix *Chandle,            // output matrix
    const GrB_Matrix M,             // optional matrix
    const GrB_Matrix A_in,          // input matrix A
    const GrB_Matrix B_in,          // input matrix B
    const GrB_Semiring semiring,    // semiring that defines C=A*B
    const bool flipxy               // if true, do z=fmult(b,a) vs fmult(a,b)
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    ASSERT_OK_OR_NULL (GB_check (M, "M for numeric C<M>=A*B", D0)) ;
    ASSERT_OK (GB_check (A_in, "A for Gustavson C=A*B", D0)) ;
    ASSERT_OK (GB_check (B_in, "B for Gustavson C=A*B", D0)) ;
    ASSERT (!PENDING (M))    ; ASSERT (!ZOMBIES (M)) ;
    ASSERT (!PENDING (A_in)) ; ASSERT (!ZOMBIES (A_in)) ;
    ASSERT (!PENDING (B_in)) ; ASSERT (!ZOMBIES (B_in)) ;
    ASSERT (A_in->vdim == B_in->vlen) ;
    ASSERT_OK (GB_check (semiring, "semiring for numeric A*B", D0)) ;

    //--------------------------------------------------------------------------
    // determine size and hypersparsity of C
    //--------------------------------------------------------------------------

    GrB_Info info ;

    (*Chandle) = NULL ;

    GrB_Matrix A2 = NULL ;
    GrB_Matrix B2 = NULL ;
    GrB_Matrix A = A_in ;
    GrB_Matrix B = B_in ;

    int64_t cvlen = A->vlen ;
    int64_t cvdim = B->vdim ;

    int64_t *Mark = NULL ;

    if (M == NULL)
    { 
        // ensure Mark is at least of size cvlen+1.  the Mark array is not used
        // for the masked matrix multiply.  This is costly if C is hypersparse,
        // but in that case the heap method will be used instead.
        OK (GB_Mark_walloc (cvlen+1)) ;  // OK: not used if C hypersparse
        Mark = GB_thread_local.Mark ;
    }

    //--------------------------------------------------------------------------
    // estimate NNZ(C) and allocate C (just the pattern)
    //--------------------------------------------------------------------------

    OK (GB_AxB_alloc (Chandle, GrB_BOOL, cvlen, cvdim, M, A, B, false,
        cvlen + NNZ (A) + NNZ (B))) ;

    GrB_Matrix C = (*Chandle) ;
    ASSERT (C != NULL) ;
    ASSERT (C->x == NULL) ;

    //==========================================================================
    // symbolic analysis when no mask is present
    //==========================================================================

    if (M == NULL)
    {
        bool A_is_hyper = IS_HYPER (A) ;
        if (A_is_hyper || IS_HYPER (B) || IS_HYPER (C))
        { 
            // symbolic analysis when one or more matrix is hypersparse
            #define HYPER
            #include "GB_AxB_Gustavson_symbolic.c"
            #undef HYPER
        }
        else
        { 
            // symbolic analysis when no matrix is hypersparse
            #include "GB_AxB_Gustavson_symbolic.c"
        }
    }

    //==========================================================================
    // numerical phase
    //==========================================================================

    // get the semiring operators
    GrB_BinaryOp mult = semiring->multiply ;
    GrB_Monoid add = semiring->add ;

    // flipxy: if true, then compute z = fmult (B(k,j), A(i,k)) instead of the
    // usual z = fmult (A(i,k), B(k,j)), since A and B have been swapped on
    // input.

    // these conditions have already been checked in the caller
    if (flipxy)
    { 
        // z = fmult (b,a) will be computed
        ASSERT (GB_Type_compatible (A->type, mult->ytype)) ;
        ASSERT (GB_Type_compatible (B->type, mult->xtype)) ;
    }
    else
    { 
        // z = fmult (a,b) will be computed
        ASSERT (GB_Type_compatible (A->type, mult->xtype)) ;
        ASSERT (GB_Type_compatible (B->type, mult->ytype)) ;
    }

    // these asserts hold for any valid semiring:
    ASSERT (mult->ztype == add->op->ztype) ;
    ASSERT (add->op->ztype == add->op->xtype) ;
    ASSERT (add->op->ztype == add->op->ytype) ;

    //--------------------------------------------------------------------------
    // allocate C->x and ensure workspace is large enough
    //--------------------------------------------------------------------------

    // C has the same type as z for z=fmult(x,y).  The type is also the
    // same as the monoid of the semiring.

    C->type = mult->ztype ;
    size_t zsize = mult->ztype->size ;

    char zwork [zsize] ;
    char cwork [zsize] ;

    GB_MALLOC_MEMORY (C->x, C->nzmax, zsize) ;
    if (C->x == NULL)
    { 
        // out of memory
        double memory = GBYTES (C->nzmax, zsize) ;
        FREE_ALL ;
        return (OUT_OF_MEMORY (memory)) ;
    }

    C->x_shallow = false ;

    // this workspace is costly if the matrices are hypersparse, but in that
    // case the heap method is used instead.
    OK (GB_Work_walloc (cvlen, zsize)) ;     // OK: not used if C hypersparse
    int8_t *Flag = NULL ;
    if (M != NULL)
    { 
        // allocate Flag
        OK (GB_Flag_walloc (cvlen)) ;        // OK: not used if C hypersparse
    }

    // w has size cvlen+1, each entry of size zsize.  Not initialized.
    void *w = GB_thread_local.Work ;

    // get the Flag array and ensure that is cleared
    if (M != NULL)
    { 
        ASSERT_FLAG_IS_CLEAR ;
        Flag = GB_thread_local.Flag ;
        ASSERT (M->vlen == C->vlen && M->vdim == C->vdim) ;
    }

    //--------------------------------------------------------------------------
    // determine the type of A2 and B2, for typecasting
    //--------------------------------------------------------------------------

    GrB_Type atype_required, btype_required ;

    if (flipxy)
    { 
        // A is passed as y, and B as x, in z = mult(x,y)
        atype_required = mult->ytype ;
        btype_required = mult->xtype ;
    }
    else
    { 
        // A is passed as x, and B as y, in z = mult(x,y)
        atype_required = mult->xtype ;
        btype_required = mult->ytype ;
    }

    //--------------------------------------------------------------------------
    // cast A and B to x and y for z=mult(x,y), if needed
    //--------------------------------------------------------------------------

    // Shallow casting creates a new matrix that has shallow pointers to the
    // prior A_in->p and A_in->i, but it constructs a new array A2->x for the
    // numerical values.  The types A_in->type and A2->type can differ, as can
    // A_in->x and A2->x, but all other content is the same.  If the types of
    // A_in and A2 are the same, then A2->x is a shallow copy of A->x and no
    // data is moved.  The CSR/CSC format of A2 and B2 is not relevant, so
    // they are kept the same as A_in and B_in.

    OK (GB_shallow_cast (&A2, atype_required, A_in->is_csc, A_in)) ;
    OK (GB_shallow_cast (&B2, btype_required, B_in->is_csc, B_in)) ;

    A = A2 ;
    B = B2 ;

    // A and B are now the right types for the multiply operator.
    // no further typecasting is needed.
    ASSERT (A->type == atype_required) ;
    ASSERT (B->type == btype_required) ;

    //--------------------------------------------------------------------------
    // compute C = A*B for built-in types and operators
    //--------------------------------------------------------------------------

    ASSERT_OK (GB_check (A->type, "A type for builtin", D0)) ;
    ASSERT_OK (GB_check (B->type, "B type for builtin", D0)) ;
    ASSERT_OK (GB_check (C->type, "C type for builtin", D0)) ;
    ASSERT_OK (GB_check (semiring, "semiring for builtin", D0)) ;

#ifndef GBCOMPACT

    // If the GB_AxB_Gustavson_builtin function has a worker for the particular
    // semiring, then it does the computation and returns done = true.
    // Otherwise, it returns done as false, and the generic worker below does
    // the work.

    // If GBCOMPACT is enabled at compile-time, then no built-in workers are
    // created, and this function is not used.  All C=A*B computations are done
    // with the generic worker below.

    bool done = false ;
    info = GB_AxB_Gustavson_builtin (C, M, A, B, semiring, flipxy, &done) ;
    ASSERT (info == GrB_SUCCESS) ;
    if (done)
    { 
        // C = A*B has been done via a hard-coded case
        ASSERT_OK (GB_check (C, "C hard-coded for numeric C=A*B", D0)) ;
        ASSERT (*Chandle == C) ;
        FREE_WORK ;
        // the masked version uses and then clears the Flag
        if (M != NULL) ASSERT_FLAG_IS_CLEAR ;
        return (REPORT_SUCCESS) ;
    }

#endif

    //--------------------------------------------------------------------------
    // generic Gustavson C=A*B for any valid semiring, built-in or user-defined
    //--------------------------------------------------------------------------

    // Define operations for GB_AxB_Gustavson_mask and GB_AxB_Gustavson_nomask

    #define IDENTITY \
        identity

    // x [i] = y
    #define COPY_SCALAR_TO_ARRAY(x,i,y,s)           \
        memcpy (x +((i)*s), y, s) ;

    // x = y [i]
    #define COPY_ARRAY_TO_SCALAR(x,y,i,s)           \
        memcpy (x, y +((i)*s), s) ;

    // x [i] = y [j]
    #define COPY_ARRAY_TO_ARRAY(x,i,y,j,s)          \
        memcpy (x +((i)*s), y +((j)*s), s);

    // generic multiply-add operation (with no mask)
    #define MULTADD_NOMASK                          \
    {                                               \
        /* w [i] += A(i,k) * B(k,j) */              \
        if (flipxy)                                 \
        {                                           \
            /* zwork = bkj * A(i,k) */              \
            fmult (zwork, bkj, Ax +(pA*asize)) ;    \
        }                                           \
        else                                        \
        {                                           \
            /* zwork = A(i,k) * bkj */              \
            fmult (zwork, Ax +(pA*asize), bkj) ;    \
        }                                           \
        /* cwork = w [i] */                         \
        memcpy (cwork, w +(i*zsize), zsize) ;       \
        /* w [i] = cwork + zwork */                 \
        fadd (w +(i*zsize), cwork, zwork) ;         \
    }

    // generic multiply-add operation (with mask)
    #define MULTADD_WITH_MASK                       \
    {                                               \
        /* w [i] += A(i,k) * B(k,j) */              \
        if (flipxy)                                 \
        {                                           \
            /* zwork = bkj * A(i,k) */              \
            fmult (zwork, bkj, Ax +(pA*asize)) ;    \
        }                                           \
        else                                        \
        {                                           \
            /* zwork = A(i,k) * bkj */              \
            fmult (zwork, Ax +(pA*asize), bkj) ;    \
        }                                           \
        if (flag > 0)                               \
        {                                           \
            /* first time C(i,j) seen */            \
            Flag [i] = -1 ;                         \
            /* w [i] = zwork */                     \
            memcpy (w +(i*zsize), zwork, zsize) ;   \
        }                                           \
        else                                        \
        {                                           \
            /* C(i,j) seen before, update it */     \
            /* cwork = w [i] */                     \
            memcpy (cwork, w +(i*zsize), zsize) ;   \
            /* w [i] = cwork + zwork */             \
            fadd (w +(i*zsize), cwork, zwork) ;     \
        }                                           \
    }

    // asize is the size of x, or y if flipxy is true, for z=mult(x,y)
    // bsize is the size of y, or x if flipxy is true, for z=mult(x,y)
    size_t asize = atype_required->size ;
    size_t bsize = btype_required->size ;

    char bkj [bsize] ;

    GB_binary_function fmult = mult->function ;
    GB_binary_function fadd  = add->op->function ;

    void *restrict identity = add->identity ;

    // get A, B, and C numerical components
    void *restrict Cx = C->x ;
    const void *restrict Ax = A->x ;
    const void *restrict Bx = B->x ;

    #include "GB_AxB_Gustavson_meta.c"

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    FREE_WORK ;     // free A2 and B2

    // cannot fail since C->plen is the upper bound: # non-empty columns of B
    ASSERT (info == GrB_SUCCESS) ;
    // if it could fail, do this:
    // OK (info) ;     // check result and return if an error occurred

    ASSERT_OK (GB_check (C, "C output for numeric C=A*B", D0)) ;
    ASSERT (*Chandle == C) ;
    return (REPORT_SUCCESS) ;
}

#undef IDENTITY
#undef COPY_SCALAR_TO_ARRAY
#undef COPY_ARRAY_TO_SCALAR
#undef COPY_ARRAY_TO_ARRAY
#undef MULTADD_NOMASK
#undef MULTADD_WITH_MASK
#undef FREE_WORK
#undef FREE_ALL
#undef OK

