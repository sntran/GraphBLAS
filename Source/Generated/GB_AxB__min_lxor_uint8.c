
//------------------------------------------------------------------------------
// GB_AxB:  hard-coded C=A*B and C<M>=A*B
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2018, All Rights Reserved.
// http://suitesparse.com   See GraphBLAS/Doc/License.txt for license.

// If this filename has a double underscore in its name ("__") then it has
// been automatically constructed from Generator/GB*AxB.[ch], via the axb*.m
// scripts, and should not be editted.  Edit the original source file instead.

//------------------------------------------------------------------------------

#include "GB.h"
#ifndef GBCOMPACT
#include "GB_heap.h"
#include "GB_AxB__semirings.h"

// The C=A*B semiring is defined by the following types and operators:

// A*B function (Gustavon):  GB_AgusB__min_lxor_uint8
// A'*B function (dot):      GB_AdotB__min_lxor_uint8
// A*B function (heap):      GB_AheapB__min_lxor_uint8
// Z type :  uint8_t (the type of C)
// XY type:  uint8_t (the type of A and B)
// Identity: UINT8_MAX (where cij = IMIN (cij,UINT8_MAX) does not change cij)
// Multiply: t = ((aik != 0) != (bkj != 0))
// Add:      cij = IMIN (cij,t)

//------------------------------------------------------------------------------
// C<M>=A*B and C=A*B: gather/scatter saxpy-based method (Gustavson)
//------------------------------------------------------------------------------

#define IDENTITY \
    UINT8_MAX

// x [i] = y
#define COPY_SCALAR_TO_ARRAY(x,i,y,s)                   \
    x [i] = y ;

// x = y [i]
#define COPY_ARRAY_TO_SCALAR(x,y,i,s)                   \
    uint8_t x = y [i] ;

// x [i] = y [i]
#define COPY_ARRAY_TO_ARRAY(x,i,y,j,s)                  \
    x [i] = y [j] ;

// multiply-add operation (no mask)
#define MULTADD_NOMASK                                  \
{                                                       \
    /* w [i] += A(i,k) * B(k,j) */                      \
    uint8_t aik = Ax [pA] ;                              \
    uint8_t t = (aik != 0) != (bkj != 0) ;                          \
    w [i] = IMIN (w [i],t) ;                                     \
}

// multiply-add operation (with mask)
#define MULTADD_WITH_MASK                               \
{                                                       \
    /* w [i] += A(i,k) * B(k,j) */                      \
    uint8_t aik = Ax [pA] ;                              \
    uint8_t t = (aik != 0) != (bkj != 0) ;                          \
    if (flag > 0)                                       \
    {                                                   \
        /* first time C(i,j) seen */                    \
        Flag [i] = -1 ;                                 \
        w [i] = t ;                                     \
    }                                                   \
    else                                                \
    {                                                   \
        /* C(i,j) seen before, update it */             \
        w [i] = IMIN (w [i],t) ;                                 \
    }                                                   \
}

GrB_Info GB_AgusB__min_lxor_uint8
(
    GrB_Matrix C,
    const GrB_Matrix M,
    const GrB_Matrix A,
    const GrB_Matrix B,
    bool flip                   // if true, A and B have been swapped
)
{ 

    uint8_t *restrict w = GB_thread_local.Work ;  // size C->vlen * zsize
    uint8_t *restrict Cx = C->x ;
    const uint8_t *restrict Ax = A->x ;
    const uint8_t *restrict Bx = B->x ;

    GrB_Info info = GrB_SUCCESS ;

    #include "GB_AxB_Gustavson_meta.c"

    return (info) ;
}

//------------------------------------------------------------------------------
// C<M>=A'*B or C=A'*B: dot product
//------------------------------------------------------------------------------

// get A(k,i)
#define DOT_GETA(pA)                                    \
    uint8_t aki = Ax [pA] ;

// get B(k,j)
#define DOT_GETB(pB)                                    \
    uint8_t bkj = Bx [pB] ;

// multiply aki*bkj, reversing them if flip is true
#define DOT_MULT(bkj)                                   \
    uint8_t t = (aki != 0) != (bkj != 0) ;

// cij += t
#define DOT_ADD                                         \
    cij = IMIN (cij,t) ;

// cij = t
#define DOT_COPY                                        \
    cij = t ;

// cij is not a pointer but a scalar; nothing to do
#define DOT_REACQUIRE ;

// clear cij
#define DOT_CLEAR                                       \
    cij = UINT8_MAX ;

// save the value of C(i,j)
#define DOT_SAVE                                        \
    Cx [cnz] = cij ;

#define DOT_WORK_TYPE \
    uint8_t

#define DOT_WORK(k) Work [k]

// Work [k] = Bx [pB]
#define DOT_SCATTER \
    Work [k] = Bx [pB] ;

GrB_Info GB_AdotB__min_lxor_uint8
(
    GrB_Matrix *Chandle,
    const GrB_Matrix M,
    const GrB_Matrix A,
    const GrB_Matrix B,
    bool flip                   // if true, A and B have been swapped
)
{ 

    GrB_Matrix C = (*Chandle) ;
    uint8_t *restrict Cx = C->x ;
    const uint8_t *restrict Ax = A->x ;
    const uint8_t *restrict Bx = B->x ;
    uint8_t cij ;
    size_t bkj_size = sizeof (uint8_t) ;

    GrB_Info info = GrB_SUCCESS ;

    #include "GB_AxB_dot_meta.c"

    return (info) ;
}

//------------------------------------------------------------------------------
// C<M>=A*B and C=A*B: heap saxpy-based method
//------------------------------------------------------------------------------

#define CIJ_GETB(pB)                                   \
    uint8_t bkj = Bx [pB] ;

// C(i,j) = A(i,k) * bkj
#define CIJ_MULT(pA)                                   \
{                                                      \
    uint8_t aik = Ax [pA] ;                             \
    cij = (aik != 0) != (bkj != 0) ;                             \
}

// C(i,j) += A(i,k) * B(k,j)
#define CIJ_MULTADD(pA,pB)                             \
{                                                      \
    uint8_t aik = Ax [pA] ;                             \
    uint8_t bkj = Bx [pB] ;                             \
    uint8_t t = (aik != 0) != (bkj != 0) ;                         \
    cij = IMIN (cij,t) ;                                      \
}

// cij is not a pointer but a scalar; nothing to do
#define CIJ_REACQUIRE ;

// cij = UINT8_MAX
#define CIJ_CLEAR                                      \
    cij = UINT8_MAX ;

// save the value of C(i,j)
#define CIJ_SAVE                                       \
    Cx [cnz] = cij ;

GrB_Info GB_AheapB__min_lxor_uint8
(
    GrB_Matrix *Chandle,
    const GrB_Matrix M,
    const GrB_Matrix A,
    const GrB_Matrix B,
    bool flip,                  // if true, A and B have been swapped
    int64_t *restrict List,
    GB_pointer_pair *restrict pA_pair,
    Element *restrict Heap,
    const int64_t bjnz_max
)
{ 

    GrB_Matrix C = (*Chandle) ;
    uint8_t *restrict Cx = C->x ;
    const uint8_t *restrict Ax = A->x ;
    const uint8_t *restrict Bx = B->x ;
    uint8_t cij ;
    int64_t cvlen = C->vlen ;

    GrB_Info info = GrB_SUCCESS ;

    #include "GB_AxB_heap_meta.c"

    return (info) ;
}

#endif

