//------------------------------------------------------------------------------
// GB_AxB_heap: compute C<M> = A*B using a heap-based method
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2018, All Rights Reserved.
// http://suitesparse.com   See GraphBLAS/Doc/License.txt for license.

//------------------------------------------------------------------------------

#include "GB.h"
#include "GB_heap.h"
#ifndef GBCOMPACT
#include "GB_AxB__semirings.h"
#endif

GrB_Info GB_AxB_heap                // C<M>=A*B or C=A*B using a heap
(
    GrB_Matrix *Chandle,            // output matrix
    const GrB_Matrix M,             // mask matrix for C<M>=A*B
    const GrB_Matrix A,             // input matrix
    const GrB_Matrix B,             // input matrix
    const GrB_Semiring semiring,    // semiring that defines C=A*B
    const bool flipxy,              // if true, do z=fmult(b,a) vs fmult(a,b)
    const int64_t bjnz_max          // max # entries in any vector of B
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    ASSERT_OK_OR_NULL (GB_check (M, "M for heap A*B", D0)) ;
    ASSERT_OK (GB_check (A, "A for heap A*B", D0)) ;
    ASSERT_OK (GB_check (B, "B for heap A*B", D0)) ;
    ASSERT (!PENDING (M)) ; ASSERT (!ZOMBIES (M)) ;
    ASSERT (!PENDING (A)) ; ASSERT (!ZOMBIES (A)) ;
    ASSERT (!PENDING (B)) ; ASSERT (!ZOMBIES (B)) ;
    ASSERT_OK (GB_Semiring_check (semiring, "semiring for heap A*B", D0)) ;
    ASSERT (A->vdim == B->vlen) ;

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
    // allocate workspace
    //--------------------------------------------------------------------------

    GrB_Info info ;

    // allocate the following arrays, all have size bjnz_max
    // (the +1 is for the mask M)
    // int64_t List [0..bjnz_max-1] ;
    // GB_pointer_pair pA_pair [0..bjnz_max1] ;
    // Element Heap [1..bjnz_max] ;

    info = GB_Work_walloc (bjnz_max,
        sizeof (int64_t) + sizeof (GB_pointer_pair) + sizeof (Element)) ;
    if (info != GrB_SUCCESS)
    { 
        // out of memory
        return (info) ;
    }

    char *w = GB_thread_local.Work ;
    int64_t *List = (int64_t *) w ;
    w += bjnz_max * sizeof (int64_t) ;
    GB_pointer_pair *pA_pair = (GB_pointer_pair *) w ;
    w += bjnz_max * sizeof (GB_pointer_pair) ;
    Element *Heap = (Element *) w ;
    Heap-- ; // Heap [0] is not used, only Heap [1..bjnz_max]

    //--------------------------------------------------------------------------
    // esimate NNZ(C) and allocate C
    //--------------------------------------------------------------------------

    int64_t cvlen = A->vlen ;
    int64_t cvdim = B->vdim ;
    GrB_Type ctype = semiring->add->op->ztype ;

    info = GB_AxB_alloc (Chandle, ctype, cvlen, cvdim, M, A, B, true,
        15 + NNZ (A) + NNZ (B)) ;

    if (info != GrB_SUCCESS)
    { 
        // out of memory
        GB_wfree ( ) ;
        return (info) ;
    }

    GrB_Matrix C = (*Chandle) ;

    //--------------------------------------------------------------------------
    // C = A*B with a heap and builtin semiring
    //--------------------------------------------------------------------------

    bool done = false ;

#ifndef GBCOMPACT

    //--------------------------------------------------------------------------
    // define the worker for the switch factory
    //--------------------------------------------------------------------------

    #define GB_AheapB(add,mult,xyname) GB_AheapB_ ## add ## mult ## xyname

    #define AxB(add,mult,xyname)                                    \
    {                                                               \
        info = GB_AheapB (add,mult,xyname) (Chandle, M, A, B,       \
            flipxy, List, pA_pair, Heap, bjnz_max) ;                \
        done = true ;                                               \
    }                                                               \
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
    // C = A*B, with a heap, and typecasting
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
        // flipxy false: aik = (xtype) A(i,k) and bkj = (ytype) B(k,j)
        // flipxy true:  aik = (ytype) A(i,k) and bkj = (xtype) B(k,j)
        char aik [flipxy ? ysize : xsize] ;
        char bkj [flipxy ? xsize : ysize] ;
        char zwork [csize] ;
        char cwork [csize] ;

        const void *Ax = A->x ;
        const void *Bx = B->x ;
        void *Cx = C->x ;
        void *cij = Cx ;        // advances through each entry of C

        void *identity = add->identity ;

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
        // C = A*B via the heap, function pointers, and typecasting
        //----------------------------------------------------------------------

        // bkj = B(k,j), located in Bx [pB]
        #define CIJ_GETB(pB)                                            \
        {                                                               \
            cast_B (bkj, Bx +((pB)*bsize), bsize) ;                     \
        }

        // C(i,j) = A(i,k) * bkj
        #define CIJ_MULT(pA)                                            \
        {                                                               \
            /* aik = A(i,k), located in Ax [pA] */                      \
            cast_A (aik, Ax +((pA)*asize), asize) ;                     \
            /* cij = aik*bkj, reversing them if flipxy is true */       \
            if (flipxy)                                                 \
            {                                                           \
                /* cij = bkj * aik */                                   \
                fmult (cij, bkj, aik) ;                                 \
            }                                                           \
            else                                                        \
            {                                                           \
                /* cij = aik * bkj */                                   \
                fmult (cij, aik, bkj) ;                                 \
            }                                                           \
        }

        // C(i,j) += A(i,k) * B(k,j)
        #define CIJ_MULTADD(pA,pB)                                      \
        {                                                               \
            /* aik = A(i,k), located in Ax [pA] */                      \
            cast_A (aik, Ax +((pA)*asize), asize) ;                     \
            /* bkj = B(k,j), located in Bx [pB] */                      \
            CIJ_GETB (pB) ;                                             \
            /* zwork = aik*bkj, reversing them if flipxy is true */     \
            if (flipxy)                                                 \
            {                                                           \
                /* zwork = bkj * aik */                                 \
                fmult (zwork, bkj, aik) ;                               \
            }                                                           \
            else                                                        \
            {                                                           \
                /* zwork = aik * bkj */                                 \
                fmult (zwork, aik, bkj) ;                               \
            }                                                           \
            /* cwork = cij */                                           \
            memcpy (cwork, cij, csize) ;                                \
            /* cij = cwork + zwork */                                   \
            fadd (cij, cwork, zwork) ;                                  \
        }

        // C->x has moved so the pointer to cij needs to be recomputed
        #define CIJ_REACQUIRE                                           \
        {                                                               \
            cij = Cx + cnz * csize ;                                    \
        }

        // cij = identity
        #define CIJ_CLEAR  memcpy (cij, identity, csize) ;

        // save the value of C(i,j) by advancing the cij pointer to next value
        #define CIJ_SAVE   cij += csize ;

        #include "GB_AxB_heap_meta.c"
    }

    //--------------------------------------------------------------------------
    // trim the size of C: this cannot fail
    //--------------------------------------------------------------------------

    info = GB_ix_realloc (C, NNZ (C), true) ;
    ASSERT (info == GrB_SUCCESS) ;
    ASSERT_OK (GB_check (C, "heap: C = A*B output", D0)) ;
    ASSERT (*Chandle == C) ;
    return (REPORT_SUCCESS) ;
}

#undef AxB
#undef GB_AheapB

#undef CIJ_GETB
#undef CIJ_MULT
#undef CIJ_MULTADD
#undef CIJ_REACQUIRE
#undef CIJ_CLEAR
#undef CIJ_SAVE

