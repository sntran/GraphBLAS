//------------------------------------------------------------------------------
// GB_AxB_dot: compute C<M> = A'*B without forming A' via dot products
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2018, All Rights Reserved.
// http://suitesparse.com   See GraphBLAS/Doc/License.txt for license.

//------------------------------------------------------------------------------

// GB_AxB_dot computes the matrix multiplication C<M>=A'*B without forming
// A' explicitly.  It is useful when A is very tall and thin (n-by-1 in
// particular).  In that case A' is costly to transpose, but A'*B is very
// easy if B is also tall and thin (say also n-by-1).

// If M is NULL, the method computes C=A'*B by considering each entry C(i,j),
// taking O(m*n) time if C is m-by-n.  This is suitable only when C is small
// (such as a scalar, a small matrix, or a vector).  If M is present, the upper
// bound on the number of entries in C is the same as nnz(M), so that space is
// allocated for C, and C(i,j) is computed only where M(i,j)=1.  This function
// assumes the mask M is not complemented.

// Compare this function with GB_AxB_Gustavson, which computes C=A*B and
// C<M>=A*B.  The computation of C=A*B requires C->p and C->i to be constructed
// first, in a symbolic phase.  Otherwise they are very similar.  The dot
// product in this algorithm is very much like the merge-add in GB_add,
// except that the merge in GB_add produces a column (a(:,j)+b(:,j)),
// whereas the merge in this function produces a scalar (a(:,j)'*b(:,j)).

#include "GB.h"
#ifndef GBCOMPACT
#include "GB_heap.h"
#include "GB_AxB__semirings.h"
#endif

GrB_Info GB_AxB_dot                 // C = A'*B using dot product method
(
    GrB_Matrix *Chandle,            // output matrix
    const GrB_Matrix M,             // mask matrix for C<M>=A'*B
    const GrB_Matrix A,             // input matrix
    const GrB_Matrix B,             // input matrix
    const GrB_Semiring semiring,    // semiring that defines C=A*B
    const bool flipxy               // if true, do z=fmult(b,a) vs fmult(a,b)
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    ASSERT_OK_OR_NULL (GB_check (M, "M for dot A'*B", D0)) ;
    ASSERT_OK (GB_check (A, "A for dot A'*B", D0)) ;
    ASSERT_OK (GB_check (B, "B for dot A'*B", D0)) ;
    ASSERT (!PENDING (M)) ; ASSERT (!ZOMBIES (M)) ;
    ASSERT (!PENDING (A)) ; ASSERT (!ZOMBIES (A)) ;
    ASSERT (!PENDING (B)) ; ASSERT (!ZOMBIES (B)) ;
    ASSERT_OK (GB_Semiring_check (semiring, "semiring for numeric A'*B", D0)) ;
    ASSERT (A->vlen == B->vlen) ;

    if (flipxy)
    { 
        // z=fmult(b,a) will be computed
        ASSERT (GB_Type_compatible (A->type, semiring->multiply->ytype)) ;
        ASSERT (GB_Type_compatible (B->type, semiring->multiply->xtype)) ;
    }
    else
    { 
        // z=fmult(a,b) will be computed
        ASSERT (GB_Type_compatible (A->type, semiring->multiply->xtype)) ;
        ASSERT (GB_Type_compatible (B->type, semiring->multiply->ytype)) ;
    }
    ASSERT (semiring->multiply->ztype == semiring->add->op->ztype) ;

    (*Chandle) = NULL ;

    //--------------------------------------------------------------------------
    // estimate NNZ(C) and allocate C
    //--------------------------------------------------------------------------

    GrB_Info info ;
    GrB_Type ctype = semiring->add->op->ztype ;
    int64_t cvlen = A->vdim ;
    int64_t cvdim = B->vdim ;

    info = GB_AxB_alloc (Chandle, ctype, cvlen, cvdim, M, A, B, true,
        15 + NNZ (A) + NNZ (B)) ;

    if (info != GrB_SUCCESS)
    { 
        // out of memory
        return (info) ;
    }

    GrB_Matrix C = (*Chandle) ;

    //--------------------------------------------------------------------------
    // C = A'*B, computing each entry with a dot product, via builtin semiring
    //--------------------------------------------------------------------------

    bool done = false ;

#ifndef GBCOMPACT

    //--------------------------------------------------------------------------
    // define the worker for the switch factory
    //--------------------------------------------------------------------------

    #define GB_AdotB(add,mult,xyname) GB_AdotB_ ## add ## mult ## xyname

    #define AxB(add,mult,xyname)                                        \
    {                                                                   \
        info = GB_AdotB (add,mult,xyname) (Chandle, M, A, B, flipxy) ;  \
        done = true ;                                                   \
    }                                                                   \
    break ;

    //--------------------------------------------------------------------------
    // launch the switch factory
    //--------------------------------------------------------------------------

    GB_Opcode mult_opcode, add_opcode ;
    GB_Type_code xycode, zcode ;

    if (GB_semiring_builtin (A, B, semiring, flipxy,
        &mult_opcode, &add_opcode, &xycode, &zcode))
    { 
        #include "GB_AxB_factory.c"
    }

    if (info != GrB_SUCCESS)
    { 
        // out of memory
        return (info) ;
    }

#endif

    //--------------------------------------------------------------------------
    // C = A'*B, computing each entry with a dot product, with typecasting
    //--------------------------------------------------------------------------

    if (!done)
    {

        //----------------------------------------------------------------------
        // get operators, functions, workspace, contents of A, B, C, and M
        //----------------------------------------------------------------------

        // get the semiring operators
        GrB_BinaryOp multiply = semiring->multiply ;
        GrB_Monoid add = semiring->add ;

        GB_binary_function fmult = multiply->function ;
        GB_binary_function fadd  = add->op->function ;

        size_t csize = C->type->size ;
        size_t asize = A->type->size ;
        size_t bsize = B->type->size ;

        size_t xsize = multiply->xtype->size ;
        size_t ysize = multiply->ytype->size ;

        // scalar workspace
        // flipxy false: aki = (xtype) A(k,i) and bkj = (ytype) B(k,j)
        // flipxy true:  aki = (ytype) A(k,i) and bkj = (xtype) B(k,j)
        size_t aki_size = flipxy ? ysize : xsize ;
        size_t bkj_size = flipxy ? xsize : ysize ;

        char aki [aki_size] ;
        char bkj [bkj_size]  ;

        char zwork [csize] ;
        char cwork [csize] ;

        const void *restrict Ax = A->x ;
        const void *restrict Bx = B->x ;
        void *restrict Cx = C->x ;
        void *restrict cij = Cx ;        // advances through each entry of C

        void *restrict identity = add->identity ;

        GB_cast_function cast_A, cast_B ;
        if (flipxy)
        { 
            // A is typecasted to y, and B is typecasted to x
            cast_A = GB_cast_factory (multiply->ytype->code, A->type->code) ;
            cast_B = GB_cast_factory (multiply->xtype->code, B->type->code) ;
        }
        else
        { 
            // A is typecasted to x, and B is typecasted to y
            cast_A = GB_cast_factory (multiply->xtype->code, A->type->code) ;
            cast_B = GB_cast_factory (multiply->ytype->code, B->type->code) ;
        }

        //----------------------------------------------------------------------
        // C = A'*B via dot products, function pointers, and typecasting
        //----------------------------------------------------------------------

        // get A(k,i)
        #define DOT_GETA(pA)                                            \
        {                                                               \
            /* aki = A(k,i), located in Ax [pA] */                      \
            cast_A (aki, Ax +((pA)*asize), asize) ;                     \
        }

        // get B(k,j)
        #define DOT_GETB(pB)                                            \
        {                                                               \
            /* bkj = B(k,j), located in Bx [pB] */                      \
            cast_B (bkj, Bx +((pB)*bsize), bsize) ;                     \
        }

        // multiply aki*bkj, reversing them if flipxy is true
        #define DOT_MULT(bkj)                                           \
        {                                                               \
            if (flipxy)                                                 \
            {                                                           \
                /* zwork = bkj * aki */                                 \
                fmult (zwork, bkj, aki) ;                               \
            }                                                           \
            else                                                        \
            {                                                           \
                /* zwork = aki * bkj */                                 \
                fmult (zwork, aki, bkj) ;                               \
            }                                                           \
        }

        // cij += zwork
        #define DOT_ADD                                                 \
        {                                                               \
            /* cwork = cij */                                           \
            memcpy (cwork, cij, csize) ;                                \
            /* cij = cwork + zwork */                                   \
            fadd (cij, cwork, zwork) ;                                  \
        }

        // cij = zwork
        #define DOT_COPY    memcpy (cij, zwork, csize) ;

        // C->x has moved so the pointer to cij needs to be recomputed
        #define DOT_REACQUIRE                                           \
        {                                                               \
            cij = Cx + cnz * csize ;                                    \
        }

        // cij = identity
        #define DOT_CLEAR   memcpy (cij, identity, csize) ;

        // save the value of C(i,j) by advancing cij pointer to next value
        #define DOT_SAVE    cij += csize ;

        #define DOT_WORK_TYPE void

        #define DOT_WORK(k) (Work +((k)*bkj_size))

        // Work [k] = (flipxy ? xtype : ytype) B (k,j)
        #define DOT_SCATTER                                             \
        {                                                               \
            /* Work [k] = B(k,j), located in Bx [p] */                  \
            cast_B (DOT_WORK (k), Bx +((pB)*bsize), bsize) ;            \
        }

        #include "GB_AxB_dot_meta.c"
    }

    //--------------------------------------------------------------------------
    // trim the size of C: this cannot fail
    //--------------------------------------------------------------------------

    info = GB_ix_realloc (C, NNZ (C), true) ;
    ASSERT (info == GrB_SUCCESS) ;
    ASSERT_OK (GB_check (C, "dot: C = A'*B output", D0)) ;
    ASSERT (*Chandle == C) ;
    return (REPORT_SUCCESS) ;
}

