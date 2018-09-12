
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

// A*B function (Gustavon):  GB_AgusB__times_min_int32
// A'*B function (dot):      GB_AdotB__times_min_int32
// A*B function (heap):      GB_AheapB__times_min_int32
// Z type :  int32_t (the type of C)
// XY type:  int32_t (the type of A and B)
// Identity: 1 (where cij *= 1 does not change cij)
// Multiply: t = (IMIN(aik,bkj))
// Add:      cij *= t

//------------------------------------------------------------------------------
// C<M>=A*B and C=A*B: gather/scatter saxpy-based method (Gustavson)
//------------------------------------------------------------------------------

#define IDENTITY \
    1

// x [i] = y
#define COPY_SCALAR_TO_ARRAY(x,i,y,s)                   \
    x [i] = y ;

// x = y [i]
#define COPY_ARRAY_TO_SCALAR(x,y,i,s)                   \
    int32_t x = y [i] ;

// x [i] = y [i]
#define COPY_ARRAY_TO_ARRAY(x,i,y,j,s)                  \
    x [i] = y [j] ;

// multiply-add operation (no mask)
#define MULTADD_NOMASK                                  \
{                                                       \
    /* w [i] += A(i,k) * B(k,j) */                      \
    int32_t aik = Ax [pA] ;                              \
    int32_t t = IMIN(aik,bkj) ;                          \
    w [i] *= t ;                                     \
}

// multiply-add operation (with mask)
#define MULTADD_WITH_MASK                               \
{                                                       \
    /* w [i] += A(i,k) * B(k,j) */                      \
    int32_t aik = Ax [pA] ;                              \
    int32_t t = IMIN(aik,bkj) ;                          \
    if (flag > 0)                                       \
    {                                                   \
        /* first time C(i,j) seen */                    \
        Flag [i] = -1 ;                                 \
        w [i] = t ;                                     \
    }                                                   \
    else                                                \
    {                                                   \
        /* C(i,j) seen before, update it */             \
        w [i] *= t ;                                 \
    }                                                   \
}

GrB_Info GB_AgusB__times_min_int32
(
    GrB_Matrix C,
    const GrB_Matrix M,
    const GrB_Matrix A,
    const GrB_Matrix B,
    bool flip                   // if true, A and B have been swapped
)
{ 

    int32_t *restrict w = GB_thread_local.Work ;  // size C->vlen * zsize
    int32_t *restrict Cx = C->x ;
    const int32_t *restrict Ax = A->x ;
    const int32_t *restrict Bx = B->x ;

    GrB_Info info = GrB_SUCCESS ;

    #include "GB_AxB_Gustavson_meta.c"

    return (info) ;
}

//------------------------------------------------------------------------------
// C<M>=A'*B or C=A'*B: dot product
//------------------------------------------------------------------------------

// get A(k,i)
#define DOT_GETA(pA)                                    \
    int32_t aki = Ax [pA] ;

// get B(k,j)
#define DOT_GETB(pB)                                    \
    int32_t bkj = Bx [pB] ;

// multiply aki*bkj, reversing them if flip is true
#define DOT_MULT(bkj)                                   \
    int32_t t = IMIN(aki,bkj) ;

// cij += t
#define DOT_ADD                                         \
    cij *= t ;

// cij = t
#define DOT_COPY                                        \
    cij = t ;

// cij is not a pointer but a scalar; nothing to do
#define DOT_REACQUIRE ;

// clear cij
#define DOT_CLEAR                                       \
    cij = 1 ;

// save the value of C(i,j)
#define DOT_SAVE                                        \
    Cx [cnz] = cij ;

#define DOT_WORK_TYPE \
    int32_t

#define DOT_WORK(k) Work [k]

// Work [k] = Bx [pB]
#define DOT_SCATTER \
    Work [k] = Bx [pB] ;

GrB_Info GB_AdotB__times_min_int32
(
    GrB_Matrix *Chandle,
    const GrB_Matrix M,
    const GrB_Matrix A,
    const GrB_Matrix B,
    bool flip                   // if true, A and B have been swapped
)
{ 

    GrB_Matrix C = (*Chandle) ;
    int32_t *restrict Cx = C->x ;
    const int32_t *restrict Ax = A->x ;
    const int32_t *restrict Bx = B->x ;
    int32_t cij ;
    size_t bkj_size = sizeof (int32_t) ;

    GrB_Info info = GrB_SUCCESS ;

    #include "GB_AxB_dot_meta.c"

    return (info) ;
}

//------------------------------------------------------------------------------
// C<M>=A*B and C=A*B: heap saxpy-based method
//------------------------------------------------------------------------------

#define CIJ_GETB(pB)                                   \
    int32_t bkj = Bx [pB] ;

// C(i,j) = A(i,k) * bkj
#define CIJ_MULT(pA)                                   \
{                                                      \
    int32_t aik = Ax [pA] ;                             \
    cij = IMIN(aik,bkj) ;                             \
}

// C(i,j) += A(i,k) * B(k,j)
#define CIJ_MULTADD(pA,pB)                             \
{                                                      \
    int32_t aik = Ax [pA] ;                             \
    int32_t bkj = Bx [pB] ;                             \
    int32_t t = IMIN(aik,bkj) ;                         \
    cij *= t ;                                      \
}

// cij is not a pointer but a scalar; nothing to do
#define CIJ_REACQUIRE ;

// cij = 1
#define CIJ_CLEAR                                      \
    cij = 1 ;

// save the value of C(i,j)
#define CIJ_SAVE                                       \
    Cx [cnz] = cij ;

GrB_Info GB_AheapB__times_min_int32
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
    int32_t *restrict Cx = C->x ;
    const int32_t *restrict Ax = A->x ;
    const int32_t *restrict Bx = B->x ;
    int32_t cij ;
    int64_t cvlen = C->vlen ;

    GrB_Info info = GrB_SUCCESS ;

    #include "GB_AxB_heap_meta.c"

    return (info) ;
}

#endif

