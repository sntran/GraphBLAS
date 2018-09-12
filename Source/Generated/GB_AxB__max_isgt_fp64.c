
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

// A*B function (Gustavon):  GB_AgusB__max_isgt_fp64
// A'*B function (dot):      GB_AdotB__max_isgt_fp64
// A*B function (heap):      GB_AheapB__max_isgt_fp64
// Z type :  double (the type of C)
// XY type:  double (the type of A and B)
// Identity: -INFINITY (where cij = FMAX (cij,-INFINITY) does not change cij)
// Multiply: t = (aik >  bkj)
// Add:      cij = FMAX (cij,t)

//------------------------------------------------------------------------------
// C<M>=A*B and C=A*B: gather/scatter saxpy-based method (Gustavson)
//------------------------------------------------------------------------------

#define IDENTITY \
    -INFINITY

// x [i] = y
#define COPY_SCALAR_TO_ARRAY(x,i,y,s)                   \
    x [i] = y ;

// x = y [i]
#define COPY_ARRAY_TO_SCALAR(x,y,i,s)                   \
    double x = y [i] ;

// x [i] = y [i]
#define COPY_ARRAY_TO_ARRAY(x,i,y,j,s)                  \
    x [i] = y [j] ;

// multiply-add operation (no mask)
#define MULTADD_NOMASK                                  \
{                                                       \
    /* w [i] += A(i,k) * B(k,j) */                      \
    double aik = Ax [pA] ;                              \
    double t = aik >  bkj ;                          \
    w [i] = FMAX (w [i],t) ;                                     \
}

// multiply-add operation (with mask)
#define MULTADD_WITH_MASK                               \
{                                                       \
    /* w [i] += A(i,k) * B(k,j) */                      \
    double aik = Ax [pA] ;                              \
    double t = aik >  bkj ;                          \
    if (flag > 0)                                       \
    {                                                   \
        /* first time C(i,j) seen */                    \
        Flag [i] = -1 ;                                 \
        w [i] = t ;                                     \
    }                                                   \
    else                                                \
    {                                                   \
        /* C(i,j) seen before, update it */             \
        w [i] = FMAX (w [i],t) ;                                 \
    }                                                   \
}

GrB_Info GB_AgusB__max_isgt_fp64
(
    GrB_Matrix C,
    const GrB_Matrix M,
    const GrB_Matrix A,
    const GrB_Matrix B,
    bool flip                   // if true, A and B have been swapped
)
{ 

    double *restrict w = GB_thread_local.Work ;  // size C->vlen * zsize
    double *restrict Cx = C->x ;
    const double *restrict Ax = A->x ;
    const double *restrict Bx = B->x ;

    GrB_Info info = GrB_SUCCESS ;

    #include "GB_AxB_Gustavson_meta.c"

    return (info) ;
}

//------------------------------------------------------------------------------
// C<M>=A'*B or C=A'*B: dot product
//------------------------------------------------------------------------------

// get A(k,i)
#define DOT_GETA(pA)                                    \
    double aki = Ax [pA] ;

// get B(k,j)
#define DOT_GETB(pB)                                    \
    double bkj = Bx [pB] ;

// multiply aki*bkj, reversing them if flip is true
#define DOT_MULT(bkj)                                   \
    double t = aki >  bkj ;

// cij += t
#define DOT_ADD                                         \
    cij = FMAX (cij,t) ;

// cij = t
#define DOT_COPY                                        \
    cij = t ;

// cij is not a pointer but a scalar; nothing to do
#define DOT_REACQUIRE ;

// clear cij
#define DOT_CLEAR                                       \
    cij = -INFINITY ;

// save the value of C(i,j)
#define DOT_SAVE                                        \
    Cx [cnz] = cij ;

#define DOT_WORK_TYPE \
    double

#define DOT_WORK(k) Work [k]

// Work [k] = Bx [pB]
#define DOT_SCATTER \
    Work [k] = Bx [pB] ;

GrB_Info GB_AdotB__max_isgt_fp64
(
    GrB_Matrix *Chandle,
    const GrB_Matrix M,
    const GrB_Matrix A,
    const GrB_Matrix B,
    bool flip                   // if true, A and B have been swapped
)
{ 

    GrB_Matrix C = (*Chandle) ;
    double *restrict Cx = C->x ;
    const double *restrict Ax = A->x ;
    const double *restrict Bx = B->x ;
    double cij ;
    size_t bkj_size = sizeof (double) ;

    GrB_Info info = GrB_SUCCESS ;

    #include "GB_AxB_dot_meta.c"

    return (info) ;
}

//------------------------------------------------------------------------------
// C<M>=A*B and C=A*B: heap saxpy-based method
//------------------------------------------------------------------------------

#define CIJ_GETB(pB)                                   \
    double bkj = Bx [pB] ;

// C(i,j) = A(i,k) * bkj
#define CIJ_MULT(pA)                                   \
{                                                      \
    double aik = Ax [pA] ;                             \
    cij = aik >  bkj ;                             \
}

// C(i,j) += A(i,k) * B(k,j)
#define CIJ_MULTADD(pA,pB)                             \
{                                                      \
    double aik = Ax [pA] ;                             \
    double bkj = Bx [pB] ;                             \
    double t = aik >  bkj ;                         \
    cij = FMAX (cij,t) ;                                      \
}

// cij is not a pointer but a scalar; nothing to do
#define CIJ_REACQUIRE ;

// cij = -INFINITY
#define CIJ_CLEAR                                      \
    cij = -INFINITY ;

// save the value of C(i,j)
#define CIJ_SAVE                                       \
    Cx [cnz] = cij ;

GrB_Info GB_AheapB__max_isgt_fp64
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
    double *restrict Cx = C->x ;
    const double *restrict Ax = A->x ;
    const double *restrict Bx = B->x ;
    double cij ;
    int64_t cvlen = C->vlen ;

    GrB_Info info = GrB_SUCCESS ;

    #include "GB_AxB_heap_meta.c"

    return (info) ;
}

#endif

