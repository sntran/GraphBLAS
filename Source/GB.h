//------------------------------------------------------------------------------
// GB.h: definitions visible only inside GraphBLAS
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2018, All Rights Reserved.
// http://suitesparse.com   See GraphBLAS/Doc/License.txt for license.

//------------------------------------------------------------------------------

// These defintions are not visible to the user.  They are used only inside
// GraphBLAS itself.

#ifndef GB_H
#define GB_H

#include "GraphBLAS.h"

//------------------------------------------------------------------------------
// code development settings
//------------------------------------------------------------------------------

// turn off debugging; do not edit these three lines
#ifndef NDEBUG
#define NDEBUG
#endif

// These flags are used for code development.  Uncomment them as needed.

// to turn on debugging, uncomment this line:
// #undef NDEBUG

// to turn on malloc debug printing, uncomment this line:
// #define PRINT_MALLOC

// to reduce code size and for faster time to compile, uncomment this line;
// GraphBLAS will be slower:
// #define GBCOMPACT

// uncomment this for code development (additional diagnostics are printed):
// #define DEVELOPER

// for coverage tests
#ifdef GBCOVER
#define GBCOVER_MAX 10000
extern int64_t gbcov [GBCOVER_MAX] ;
extern int gbcover_max ;
#endif

//------------------------------------------------------------------------------
// opaque content of GraphBLAS objects
//------------------------------------------------------------------------------

#define GB_LEN 128

struct GB_Type_opaque       // content of GrB_Type
{
    int64_t magic ;         // for detecting uninitialized objects
    size_t size ;           // size of the type
    int code ;              // the type code
    char name [GB_LEN] ;    // name of the type
} ;

struct GB_UnaryOp_opaque    // content of GrB_UnaryOp
{
    int64_t magic ;         // for detecting uninitialized objects
    GrB_Type xtype ;        // type of x
    GrB_Type ztype ;        // type of z
    void *function ;        // a pointer to the unary function
    char name [GB_LEN] ;    // name of the unary operator
    int opcode ;            // operator opcode
} ;

struct GB_BinaryOp_opaque   // content of GrB_BinaryOp
{
    int64_t magic ;         // for detecting uninitialized objects
    GrB_Type xtype ;        // type of x
    GrB_Type ytype ;        // type of y
    GrB_Type ztype ;        // type of z
    void *function ;        // a pointer to the binary function
    char name [GB_LEN] ;    // name of the binary operator
    int opcode ;            // operator opcode
} ;

struct GB_SelectOp_opaque   // content of GxB_SelectOp
{
    int64_t magic ;         // for detecting uninitialized objects
    GrB_Type xtype ;        // type of x, or NULL if generic
    void *function ;        // a pointer to the select function
    char name [GB_LEN] ;    // name of the select operator
    int opcode ;            // operator opcode
} ;

struct GB_Monoid_opaque     // content of GrB_Monoid
{
    int64_t magic ;         // for detecting uninitialized objects
    GrB_BinaryOp op ;       // binary operator of the monoid
    void *identity ;        // identity of the monoid; size is op->ztype->size
    bool user_defined ;     // true if monoid is user-defined
} ;

struct GB_Semiring_opaque   // content of GrB_Semiring
{
    int64_t magic ;         // for detecting uninitialized objects
    GrB_Monoid add ;        // add operator of the semiring
    GrB_BinaryOp multiply ; // multiply operator of the semiring
    bool user_defined ;     // true if semiring is user-defined
} ;

struct GB_Vector_opaque     // content of GrB_Vector: m-by-1 standard CSC matrix
{
    #include "GB_matrix.h"
} ;

struct GB_Matrix_opaque     // content of GrB_Matrix
{
    #include "GB_matrix.h"
} ;

struct GB_Descriptor_opaque // content of GrB_Descriptor
{
    int64_t magic ;         // for detecting uninitialized objects
    GrB_Desc_Value out ;    // output descriptor
    GrB_Desc_Value mask ;   // mask descriptor
    GrB_Desc_Value in0 ;    // first input descriptor (A for C=A*B, for example)
    GrB_Desc_Value in1 ;    // second input descriptor (B for C=A*B)
    GrB_Desc_Value axb ;    // for selecting the method for C=A*B
} ;

//------------------------------------------------------------------------------
// default options
//------------------------------------------------------------------------------

// These parameters define the content of extern const values that can be
// used as inputs to GxB_*Option_set.

// The default format is by column (CSC), with a hyper_ratio of 1/16.

#define GB_HYPER_DEFAULT (0.0625)

// compile SuiteSparse:GraphBLAS with "-DBYROW" to make GxB_BY_ROW the default
// format
#ifdef BYROW
#define GB_FORMAT_DEFAULT GxB_BY_ROW
#else
#define GB_FORMAT_DEFAULT GxB_BY_COL
#endif

// these parameters define the hyper_ratio needed to ensure matrix stays
// either always hypersparse, or never hypersparse.
#define GB_ALWAYS_HYPER (1.0)
#define GB_NEVER_HYPER  (-1.0)

#define GB_FORCE_HYPER 1
#define GB_FORCE_NONHYPER 0
#define GB_AUTO_HYPER (-1)

#define SAME_HYPER_AS(A_is_hyper) \
    ((A_is_hyper) ? GB_FORCE_HYPER : GB_FORCE_NONHYPER)

// if A is hypersparse but all vectors are present, then
// treat A as if it were non-hypersparse
#define IS_HYPER(A) \
    (((A) != NULL) && ((A)->is_hyper && ((A)->nvec < (A)->vdim)))

//------------------------------------------------------------------------------
// macros for matrices and vectors
//------------------------------------------------------------------------------

// If A->nzmax is zero, then A->p might not be allocated.
// Note that this function does not count pending tuples;
// use WAIT(A) first, if needed.
#define NNZ(A) (((A)->nzmax > 0) ? (A)->p [(A)->nvec] : 0)

// A is nrows-by-ncols, in either CSR or CSC format
#define NROWS(A) (A->is_csc ? A->vlen : A->vdim)
#define NCOLS(A) (A->is_csc ? A->vdim : A->vlen)

// The internal content of a GrB_Matrix and GrB_Vector are identical, and
// inside SuiteSparse:GraphBLAS, they can be typecasted between each other.
// This typecasting feature should not be done in user code, however, since it
// is not supported in the API.  All GrB_Vector objects can be safely
// typecasted into a GrB_Matrix, but not the other way around.  The GrB_Vector
// object is more restrictive.  The VECTOR_OK(v) macro defines the content that
// all GrB_Vector objects must have.

// VECTOR_OK(v) is used mainly for assertions, but also to determine when it is
// safe to typecast an n-by-1 GrB_Matrix (in standard CSC format) into a
// GrB_Vector.  This is not done in the main SuiteSparse:GraphBLAS library, but
// in the GraphBLAS/Test directory only.  The macro is also used in
// GB_Vector_check, to ensure the content of a GrB_Vector is valid.

#define VECTOR_OK(v)                \
(                                   \
    ((v) != NULL) &&                \
    ((v)->is_hyper == false) &&     \
    ((v)->is_csc == true) &&        \
    ((v)->plen == 1) &&             \
    ((v)->vdim == 1) &&             \
    ((v)->nvec == 1) &&             \
    ((v)->h == NULL) &&             \
    ((v)->j_pending == NULL)        \
)                                   \

//------------------------------------------------------------------------------
// GB_INDEX_MAX
//------------------------------------------------------------------------------

// The largest valid dimension permitted in this implementation is 2^60.
// Matrices with that many rows and/or columns can be actually be easily
// created, particulary if they are hypersparse since in that case O(nrows) or
// O(ncols) memory is not needed.  For the standard formats, O(ncols) space is
// needed for CSC and O(nrows) space is needed for CSR.

// The time complexity of many operations does not depend nrows or ncols at
// all.  Even some forms of matrix multiply can be performed: C=A'*B and
// w=u'*A, for example, without a time or memory complexity depending on nrows,
// for the CSC format.

// MATLAB has a limit of 2^48-1
#ifdef MATLAB_MEX_FILE
#define GB_INDEX_MAX ((GrB_Index) ((1ULL << 48)-1))
#else
#define GB_INDEX_MAX ((GrB_Index) (1ULL << 60))
#endif

// format strings, normally %llu and %lld, for GrB_Index values
#define GBu "%" PRIu64
#define GBd "%" PRId64

//------------------------------------------------------------------------------
// GraphBLAS check functions: check and optionally print an object
//------------------------------------------------------------------------------

// pr values for *_check functions
#define D3 3
#define D2 2
#define D1 1
#define D0 0

GrB_Info GB_Type_check      // check a GraphBLAS Type
(
    const GrB_Type type,    // GraphBLAS type to print and check
    const char *name,       // name of the type
    int pr                  // 0: print nothing, 1: print header and errors,
                            // 2: print brief, 3: print all
) ;

GrB_Info GB_BinaryOp_check  // check a GraphBLAS binary operator
(
    const GrB_BinaryOp op,  // GraphBLAS operator to print and check
    const char *name,       // name of the operator
    int pr                  // 0: print nothing, 1: print header and errors,
                            // 2: print brief, 3: print all
) ;

GrB_Info GB_UnaryOp_check   // check a GraphBLAS unary operator
(
    const GrB_UnaryOp op,   // GraphBLAS operator to print and check
    const char *name,       // name of the operator
    int pr                  // 0: print nothing, 1: print header and errors,
                            // 2: print brief, 3: print all
) ;

GrB_Info GB_SelectOp_check  // check a GraphBLAS select operator
(
    const GxB_SelectOp op,  // GraphBLAS operator to print and check
    const char *name,       // name of the operator
    int pr                  // 0: print nothing, 1: print header and errors,
                            // 2: print brief, 3: print all
) ;

GrB_Info GB_Monoid_check        // check a GraphBLAS monoid
(
    const GrB_Monoid monoid,    // GraphBLAS monoid to print and check
    const char *name,           // name of the monoid
    int pr                      // 0: print nothing, 1: print header and errors,
                                // 2: print brief, 3: print all
) ;

GrB_Info GB_Semiring_check          // check a GraphBLAS semiring
(
    const GrB_Semiring semiring,    // GraphBLAS semiring to print and check
    const char *name,               // name of the semiring
    int pr                          // 0: print nothing, 1: print header and
                                    // errors, 2: print brief, 3: print all
) ;

GrB_Info GB_Descriptor_check    // check a GraphBLAS descriptor
(
    const GrB_Descriptor D,     // GraphBLAS descriptor to print and check
    const char *name,           // name of the descriptor
    int pr                      // 0: print nothing, 1: print header and
                                // errors, 2: print brief, 3: print all
) ;

GrB_Info GB_matvec_check    // check a GraphBLAS matrix or vector
(
    const GrB_Matrix A,     // GraphBLAS matrix to print and check
    const char *name,       // name of the matrix
    int pr,                 // 0: print nothing, 1: print header and errors,
                            // 2: print brief, 3: print all
    const char *kind        // "matrix" or "vector"
) ;

GrB_Info GB_Matrix_check    // check a GraphBLAS matrix
(
    const GrB_Matrix A,     // GraphBLAS matrix to print and check
    const char *name,       // name of the matrix
    int pr                  // 0: print nothing, 1: print header and errors,
                            // 2: print brief, 3: print all
) ;

GrB_Info GB_Vector_check    // check a GraphBLAS vector
(
    const GrB_Vector A,     // GraphBLAS vector to print and check
    const char *name,       // name of the vector
    int pr                  // 0: print nothing, 1: print header and errors,
                            // 2: print brief, 3: print all
) ;

#define GB_check(x,name,pr)                             \
    _Generic                                            \
    (                                                   \
        (x),                                            \
        const GrB_Type       : GB_Type_check       ,    \
              GrB_Type       : GB_Type_check       ,    \
        const GrB_BinaryOp   : GB_BinaryOp_check   ,    \
              GrB_BinaryOp   : GB_BinaryOp_check   ,    \
        const GxB_SelectOp   : GB_SelectOp_check   ,    \
              GxB_SelectOp   : GB_SelectOp_check   ,    \
        const GrB_UnaryOp    : GB_UnaryOp_check    ,    \
              GrB_UnaryOp    : GB_UnaryOp_check    ,    \
        const GrB_Monoid     : GB_Monoid_check     ,    \
              GrB_Monoid     : GB_Monoid_check     ,    \
        const GrB_Semiring   : GB_Semiring_check   ,    \
              GrB_Semiring   : GB_Semiring_check   ,    \
        const GrB_Matrix     : GB_Matrix_check     ,    \
              GrB_Matrix     : GB_Matrix_check     ,    \
        const GrB_Vector     : GB_Vector_check     ,    \
              GrB_Vector     : GB_Vector_check     ,    \
        const GrB_Descriptor : GB_Descriptor_check ,    \
              GrB_Descriptor : GB_Descriptor_check      \
    ) (x, name, pr)

//------------------------------------------------------------------------------
// debugging definitions
//------------------------------------------------------------------------------

#ifdef MATLAB_MEX_FILE
// compiling GraphBLAS in a MATLAB mexFunction.  Use mxMalloc, mxFree, etc.
#include "mex.h"
#include "matrix.h"
#undef MALLOC
#undef FREE
#undef CALLOC
#undef REALLOC
#define MALLOC  mxMalloc
#define FREE    mxFree
#define CALLOC  mxCalloc
#define REALLOC mxRealloc
#define malloc  mxMalloc
#define free    mxFree
#define calloc  mxCalloc
#define realloc mxRealloc
#endif

#ifndef NDEBUG

    // debugging enabled
    #ifdef MATLAB_MEX_FILE
    #define ASSERT(x) \
    {                                                                       \
        if (!(x))                                                           \
        {                                                                   \
            mexErrMsgTxt ("failure: " __FILE__ " line: " GB_XSTR(__LINE__)) ; \
        }                                                                   \
    }
    #else
    #include <assert.h>
    #define ASSERT(x) assert (x) ;
    #endif
    #define ASSERT_OK(X)                                                    \
    {                                                                       \
        GrB_Info Info = X ;                                                 \
        ASSERT (Info == GrB_SUCCESS) ;                                      \
    }
    #define ASSERT_OK_OR_NULL(X)                                            \
    {                                                                       \
        GrB_Info Info = X ;                                                 \
        ASSERT (Info == GrB_SUCCESS || Info == GrB_NULL_POINTER) ;          \
    }
    #define ASSERT_OK_OR_JUMBLED(X)                                         \
    {                                                                       \
        GrB_Info Info = X ;                                                 \
        ASSERT (Info == GrB_SUCCESS || Info == GrB_INDEX_OUT_OF_BOUNDS) ;   \
    }

#else

    // debugging disabled
    #define ASSERT(x)
    #define ASSERT_OK(X)
    #define ASSERT_OK_OR_NULL(X)
    #define ASSERT_OK_OR_JUMBLED(X)

#endif

#define IMPLIES(p,q) (!(p) || (q))

// for finding tests that trigger statement coverage
#ifdef MATLAB_MEX_FILE
#define GOTCHA \
mexErrMsgTxt ("gotcha: " __FILE__ " line: " GB_XSTR(__LINE__)) ;
#else
#define GOTCHA \
{ printf ("gotcha: " __FILE__ " line: " GB_XSTR(__LINE__)"\n") ; abort () ; }
#endif

#define Here printf (" Here: %s line: %d\n",  __FILE__, __LINE__) ;

//------------------------------------------------------------------------------
// aliased objects
//------------------------------------------------------------------------------

// GraphBLAS allows all inputs to all user-accessible objects to be aliased, as
// in GrB_mxm (C, C, C, C, ...), which is valid.  Internal routines are more
// restrictive.

// true if C and X are aliased and not NULL
#define ALIASED(C,X)            ((C) == (X) && (C) != NULL)

// C may not be aliased with X, Y, or Z.  But X, Y and/or Z may be aliased
// with each other.
#define NOT_ALIASED(C,X)        (!ALIASED (C,X))
#define NOT_ALIASED_2(C,X,Y)     (NOT_ALIASED (C,X) && NOT_ALIASED  (C,Y))
#define NOT_ALIASED_3(C,X,Y,Z)   (NOT_ALIASED (C,X) && NOT_ALIASED_2 (C,Y,Z))

// these macros are always true but are used just for commenting the code
#define ALIAS_OK(C,X)           (ALIASED (C,X) || NOT_ALIASED (C,X))
#define ALIAS_OK2(C,X,Y)        (ALIAS_OK (C,X) && ALIAS_OK  (C,Y))
#define ALIAS_OK3(C,X,Y,Z)      (ALIAS_OK (C,X) && ALIAS_OK2 (C,Y,Z))

//------------------------------------------------------------------------------
// GraphBLAS memory manager
//------------------------------------------------------------------------------

// GraphBLAS can be compiled with -DMALLOC=mymallocfunc to redefine
// the malloc function and other memory management functions.
// By default, these are simply the system malloc, free, etc, routines.
#ifndef MALLOC
#define MALLOC malloc
#endif

#ifndef CALLOC
#define CALLOC calloc
#endif

#ifndef REALLOC
#define REALLOC realloc
#endif

#ifndef FREE
#define FREE free
#endif

#define GBYTES(n,s)  ((((double) (n)) * ((double) (s))) / 1e9)

//------------------------------------------------------------------------------
// internal GraphBLAS type and operator codes
//------------------------------------------------------------------------------

// MAGIC is an arbitrary number that is placed inside each object when it is
// initialized, as a way of detecting uninitialized objects.
#define MAGIC  0x00981B0787374E72

// The magic number is set to FREED when the object is freed, as a way of
// helping to detect dangling pointers.
#define FREED  0x0911911911911911

// The value is set to MAGIC2 when the object has been allocated but cannot yet
// be used in most methods and operations.  Currently this is used only for
// when A->p array is allocated but not initialized.
#define MAGIC2 0x10981B0787374E72

typedef enum
{
    // the 12 scalar types
    GB_BOOL_code,               // 0: bool
    GB_INT8_code,               // 1: int8_t
    GB_UINT8_code,              // 2: uint8_t
    GB_INT16_code,              // 3: int16_t
    GB_UINT16_code,             // 4: uint16_t
    GB_INT32_code,              // 5: int32_t
    GB_UINT32_code,             // 6: uint32_t
    GB_INT64_code,              // 7: int64_t
    GB_UINT64_code,             // 8: uint64_t
    GB_FP32_code,               // 9: float
    GB_FP64_code,               // 10: double
    GB_UDT_code,                // 11: void *, user-defined type
}
GB_Type_code ;                  // enumerated type code

// predefined type objects
extern struct GB_Type_opaque
    GB_opaque_BOOL   ,  // GrB_BOOL is a pointer to this object, etc.
    GB_opaque_INT8   ,
    GB_opaque_UINT8  ,
    GB_opaque_INT16  ,
    GB_opaque_UINT16 ,
    GB_opaque_INT32  ,
    GB_opaque_UINT32 ,
    GB_opaque_INT64  ,
    GB_opaque_UINT64 ,
    GB_opaque_FP32   ,
    GB_opaque_FP64   ;

// operator codes used in GrB_BinaryOp and GrB_UnaryOp structures
typedef enum
{
    //--------------------------------------------------------------------------
    // NOP
    //--------------------------------------------------------------------------

    // a placeholder; not an actual operator
    GB_NOP_opcode,      //  0: nop

    //--------------------------------------------------------------------------
    // T -> T
    //--------------------------------------------------------------------------

    // 6 unary operators x=f(x) that return the same type as their input
    GB_ONE_opcode,      //  1: z = 1
    GB_IDENTITY_opcode, //  2: z = x
    GB_AINV_opcode,     //  3: z = -x
    GB_ABS_opcode,      //  4: z = abs(x)
    GB_MINV_opcode,     //  5: z = 1/x ; special cases for bool and integers
    GB_LNOT_opcode,     //  6: z = !x

    //--------------------------------------------------------------------------
    // TxT -> T
    //--------------------------------------------------------------------------

    // 8 binary operators z=f(x,y) that return the same type as their inputs
    GB_FIRST_opcode,    //  7: z = x
    GB_SECOND_opcode,   //  8: z = y
    GB_MIN_opcode,      //  9: z = min(x,y)
    GB_MAX_opcode,      // 10: z = max(x,y)
    GB_PLUS_opcode,     // 11: z = x + y
    GB_MINUS_opcode,    // 12: z = x - y
    GB_TIMES_opcode,    // 13: z = x * y
    GB_DIV_opcode,      // 14: z = x / y ; special cases for bool and ints

    // 6 binary operators z=f(x,y), x,y,z all the same type
    GB_ISEQ_opcode,     // 15: z = (x == y)
    GB_ISNE_opcode,     // 16: z = (x != y)
    GB_ISGT_opcode,     // 17: z = (x >  y)
    GB_ISLT_opcode,     // 18: z = (x <  y)
    GB_ISGE_opcode,     // 19: z = (x >= y)
    GB_ISLE_opcode,     // 20: z = (x <= y)

    // 3 binary operators that work on purely boolean values
    GB_LOR_opcode,      // 21: z = (x != 0) || (y != 0)
    GB_LAND_opcode,     // 22: z = (x != 0) && (y != 0)
    GB_LXOR_opcode,     // 23: z = (x != 0) != (y != 0)

    //--------------------------------------------------------------------------
    // TxT -> bool
    //--------------------------------------------------------------------------

    // 6 binary operators z=f(x,y) that return bool (TxT -> bool)
    GB_EQ_opcode,       // 24: z = (x == y)
    GB_NE_opcode,       // 25: z = (x != y)
    GB_GT_opcode,       // 26: z = (x >  y)
    GB_LT_opcode,       // 27: z = (x <  y)
    GB_GE_opcode,       // 28: z = (x >= y)
    GB_LE_opcode,       // 29: z = (x <= y)

    //--------------------------------------------------------------------------
    // user-defined, both unary and binary
    //--------------------------------------------------------------------------

    // all user-defined operators are given this code (both unary and binary)
    GB_USER_opcode      // 30: user defined operator
}
GB_Opcode ;

#define NAME ((name != NULL) ? name : "")

#define TYPE            bool
#define GB(x)           GB_ ## x ## _BOOL
#define CAST_NAME(x)    GB_cast_bool_ ## x
#include "GB_ops_template.h"

#define TYPE            int8_t
#define GB(x)           GB_ ## x ## _INT8
#define CAST_NAME(x)    GB_cast_int8_t_ ## x
#include "GB_ops_template.h"

#define TYPE            uint8_t
#define GB(x)           GB_ ## x ## _UINT8
#define CAST_NAME(x)    GB_cast_uint8_t_ ## x
#include "GB_ops_template.h"

#define TYPE            int16_t
#define GB(x)           GB_ ## x ## _INT16
#define CAST_NAME(x)    GB_cast_int16_t_ ## x
#include "GB_ops_template.h"

#define TYPE            uint16_t
#define GB(x)           GB_ ## x ## _UINT16
#define CAST_NAME(x)    GB_cast_uint16_t_ ## x
#include "GB_ops_template.h"

#define TYPE            int32_t
#define GB(x)           GB_ ## x ## _INT32
#define CAST_NAME(x)    GB_cast_int32_t_ ## x
#include "GB_ops_template.h"

#define TYPE            uint32_t
#define GB(x)           GB_ ## x ## _UINT32
#define CAST_NAME(x)    GB_cast_uint32_t_ ## x
#include "GB_ops_template.h"

#define TYPE            int64_t
#define GB(x)           GB_ ## x ## _INT64
#define CAST_NAME(x)    GB_cast_int64_t_ ## x
#include "GB_ops_template.h"

#define TYPE            uint64_t
#define GB(x)           GB_ ## x ## _UINT64
#define CAST_NAME(x)    GB_cast_uint64_t_ ## x
#include "GB_ops_template.h"

#define TYPE            float
#define GB(x)           GB_ ## x ## _FP32
#define CAST_NAME(x)    GB_cast_float_ ## x
#include "GB_ops_template.h"

#define TYPE            double
#define GB(x)           GB_ ## x ## _FP64
#define CAST_NAME(x)    GB_cast_double_ ## x
#include "GB_ops_template.h"

//------------------------------------------------------------------------------
// select opcodes
//------------------------------------------------------------------------------

// operator codes used in GrB_SelectOp structures
typedef enum
{
    // built-in select operators:
    GB_TRIL_opcode,
    GB_TRIU_opcode,
    GB_DIAG_opcode,
    GB_OFFDIAG_opcode,
    GB_NONZERO_opcode,

    // for all user-defined select operators:
    GB_USER_SELECT_opcode
}
GB_Select_Opcode ;

//------------------------------------------------------------------------------
// internal GraphBLAS functions
//------------------------------------------------------------------------------

typedef enum                    // input parameter to GB_new and GB_create
{
    GB_Ap_calloc,               // 0: calloc A->p, malloc A->h if hypersparse
    GB_Ap_malloc,               // 1: malloc A->p, malloc A->h if hypersparse
    GB_Ap_null                  // 2: do not allocate A->p or A->h
}
GB_Ap_code ;

GrB_Info GB_new                 // create matrix, except for indices & values
(
    GrB_Matrix *matrix_handle,  // handle of matrix to create
    const GrB_Type type,        // matrix type
    const int64_t vlen,         // length of each vector
    const int64_t vdim,         // number of vectors
    const GB_Ap_code Ap_option, // calloc/malloc A->p and A->h, or leave NULL
    const bool is_csc,          // true if CSC, false if CSR
    const int hyper_option,     // 1:hyper, 0:nonhyper, -1:auto
    const double hyper_ratio,   // A->hyper_ratio, unless auto
    const int64_t plen          // size of A->p and A->h, if hypersparse
) ;

GrB_Info GB_create              // create a new matrix, including A->i and A->x
(
    GrB_Matrix *matrix_handle,  // output matrix to create
    const GrB_Type type,        // type of output matrix
    const int64_t vlen,         // length of each vector
    const int64_t vdim,         // number of vectors
    const GB_Ap_code Ap_option, // calloc/malloc A->p and A->h, or leave NULL
    const bool is_csc,          // true if CSC, false if CSR
    const int hyper_option,     // 1:hyper, 0:nonhyper, -1:auto
    const double hyper_ratio,   // A->hyper_ratio, unless auto
    const int64_t plen,         // size of A->p and A->h, if hypersparse
    const int64_t anz,          // number of nonzeros the matrix must hold
    const bool numeric          // if true, allocate A->x, else A->x is NULL
) ;

GrB_Info GB_hyper_realloc
(
    GrB_Matrix A,               // matrix with hyperlist to reallocate
    int64_t plen_new            // new size of A->p and A->h
) ;

GrB_Info GB_clear           // clear a matrix, type and dimensions unchanged
(
    GrB_Matrix A            // matrix to clear
) ;

GrB_Info GB_dup             // make an exact copy of a matrix
(
    GrB_Matrix *handle,     // handle of output matrix to create
    const GrB_Matrix A      // input matrix to copy
) ;

GrB_Info GB_nvals           // get the number of entries in a matrix
(
    GrB_Index *nvals,       // matrix has nvals entries
    const GrB_Matrix A      // matrix to query
) ;

GrB_Info GB_type            // get the type of a matrix
(
    GrB_Type *type_handle,  // returns the type of the matrix
    const GrB_Matrix A      // matrix to query
) ;

GrB_Info GB_ix_alloc        // allocate A->i and A->x space in a matrix
(
    GrB_Matrix A,           // matrix to allocate space for
    const GrB_Index nzmax,  // number of entries the matrix can hold
    const bool numeric      // if true, allocate A->x, otherwise A->x is NULL
) ;

GrB_Info GB_ix_realloc      // reallocate space in a matrix
(
    GrB_Matrix A,           // object to allocate space for
    const GrB_Index nzmax,  // new number of entries the object can hold
    const bool numeric      // if true, reallocate A->x, otherwise A->x is NULL
) ;

GrB_Info GB_ix_resize           // resize a matrix
(
    GrB_Matrix A,
    const int64_t anz_new       // required new NNZ(A)
) ;

void GB_ix_free                 // free A->i and A->x of a matrix
(
    GrB_Matrix A                // matrix with content to free
) ;

void GB_ph_free                 // free A->p and A->h of a matrix
(
    GrB_Matrix A                // matrix with content to free
) ;

void GB_phix_free               // free all content of a matrix
(
    GrB_Matrix A                // matrix with content to free
) ;

void GB_free                    // free a matrix
(
    GrB_Matrix *matrix_handle   // handle of matrix to free
) ;

bool GB_Type_compatible         // check if two types can be typecast
(
    const GrB_Type atype,
    const GrB_Type btype
) ;

bool GB_code_compatible         // check if two types can be typecast
(
    const GB_Type_code acode,
    const GB_Type_code bcode
) ;

GrB_Info GB_transpose_bucket    // bucket transpose; typecast and apply op
(
    GrB_Matrix *Chandle,        // output matrix (unallocated on input)
    const GrB_Type ctype,       // type of output matrix C
    const bool C_is_csc,        // format of output matrix C
    const GrB_Matrix A,         // input matrix
    const GrB_UnaryOp op        // operator to apply, NULL if no operator
) ;

void GB_transpose_ix        // transpose the pattern and values of a matrix
(
    int64_t *Rp,            // size m+1, input: row pointers, shifted on output
    int64_t *Ri,            // size cnz, output column indices
    void *Rx,               // size cnz, output numerical values, type R_type
    const GrB_Type R_type,  // type of output R (do typecasting into R)
    const GrB_Matrix A      // input matrix
) ;

void GB_transpose_op        // transpose and apply an operator to a matrix
(
    int64_t *Rp,            // size m+1, input: row pointers, shifted on output
    int64_t *Ri,            // size cnz, output column indices
    void *Rx,               // size cnz, output values, type op->ztype
    const GrB_UnaryOp op,   // operator to apply, NULL if no operator
    const GrB_Matrix A      // input matrix
) ;

GrB_Info GB_apply                   // C<M> = accum (C, op(A)) or op(A')
(
    GrB_Matrix C,                   // input/output matrix for results
    const bool C_replace,           // C descriptor
    const GrB_Matrix M,             // optional mask for C, unused if NULL
    const bool Mask_comp,           // M descriptor
    const GrB_BinaryOp accum,       // optional accum for Z=accum(C,T)
    const GrB_UnaryOp op,           // operator to apply to the entries
    const GrB_Matrix A,             // first input:  matrix A
    bool A_transpose                // A matrix descriptor
) ;

GrB_Info GB_select          // C<Mask> = accum (C, select(A,k)) or select(A',k)
(
    GrB_Matrix C,                   // input/output matrix for results
    const bool C_replace,           // C descriptor
    const GrB_Matrix Mask,          // optional mask for C, unused if NULL
    const bool Mask_comp,           // Mask descriptor
    const GrB_BinaryOp accum,       // optional accum for Z=accum(C,T)
    const GxB_SelectOp op,          // operator to select the entries
    const GrB_Matrix A,             // input matrix
    const void *k,                  // optional input for select operator
    bool A_transpose                // A matrix descriptor
) ;

GrB_Info GB_shallow_cast    // create a shallow typecasted matrix
(
    GrB_Matrix *Chandle,    // output matrix C, of type op->ztype
    const GrB_Type ctype,   // type of the output matrix C
    const bool C_is_csc,    // desired CSR/CSC format of C
    const GrB_Matrix A      // input matrix to typecast
) ;

GrB_Info GB_shallow_op      // create shallow matrix and apply operator
(
    GrB_Matrix *Chandle,    // output matrix C, of type op->ztype
    const bool C_is_csc,    // desired CSR/CSC format of C
    const GrB_UnaryOp op,   // operator to apply
    const GrB_Matrix A      // input matrix to typecast
) ;

void GB_cast_array              // typecast an array
(
    void *C,                    // output array
    const GB_Type_code code1,   // type code for C
    const void *A,              // input array
    const GB_Type_code code2,   // type code for A
    const int64_t n             // number of entries in C and A
) ;

typedef void (*GB_binary_function) (void *, const void *, const void *) ;

typedef void (*GB_unary_function)  (void *, const void *) ;

typedef bool (*GB_select_function)      // return true if A(i,j) is kept
(
    GrB_Index i,                // row index of A(i,j)
    GrB_Index j,                // column index of A(i,j)
    GrB_Index nrows,            // number of rows of A
    GrB_Index ncols,            // number of columns of A
    const void *x,              // value of A(i,j)
    const void *k               // optional input for select function
) ;

typedef void (*GB_cast_function)   (void *, const void *, size_t) ;

GB_cast_function GB_cast_factory   // returns pointer to function to cast x to z
(
    const GB_Type_code code1,      // the type of z, the output value
    const GB_Type_code code2       // the type of x, the input value
) ;

void GB_copy_user_user (void *z, const void *x, size_t size) ;

GrB_Info GB_add             // C = A+B
(
    GrB_Matrix *Chandle,    // output matrix (unallocated on input)
    const GrB_Type ctype,   // type of output matrix C
    const bool C_is_csc,    // format of output matrix C
    const GrB_Matrix A,     // input A matrix
    const GrB_Matrix B,     // input B matrix
    const GrB_BinaryOp op   // op to perform C = op (A,B)
) ;

GrB_Info GB_emult           // C = A.*B
(
    GrB_Matrix *Chandle,    // output matrix (unallocated on input)
    const GrB_Type ctype,   // type of output matrix C
    const bool C_is_csc,    // format of output matrix C
    const GrB_Matrix A,     // input A matrix
    const GrB_Matrix B,     // input B matrix
    const GrB_BinaryOp op   // op to perform C = op (A,B)
) ;

GrB_Info GB_transplant          // transplant one matrix into another
(
    GrB_Matrix C,               // matrix to overwrite with A
    const GrB_Type ctype,       // new type of C
    GrB_Matrix *Ahandle         // matrix to copy from and free
) ;

GrB_Info GB_transplant_conform      // transplant and conform hypersparsity
(
    GrB_Matrix C,                   // destination matrix to transplant into
    GrB_Type ctype,                 // type to cast into
    GrB_Matrix *Thandle             // source matrix
) ;

size_t GB_code_size             // return the size of a type, given its code
(
    const GB_Type_code code,    // input code of the type to find the size of
    const size_t user_size      // known size of user-defined type
) ;

void *GB_calloc_memory      // pointer to allocated block of memory
(
    size_t nitems,          // number of items to allocate
    size_t size_of_item     // sizeof each item
) ;

void *GB_malloc_memory      // pointer to allocated block of memory
(
    size_t nitems,          // number of items to allocate
    size_t size_of_item     // sizeof each item
) ;

void *GB_realloc_memory     // pointer to reallocated block of memory, or
                            // to original block if the realloc failed.
(
    size_t nitems_new,      // new number of items in the object
    size_t nitems_old,      // old number of items in the object
    size_t size_of_item,    // sizeof each item
    void *p,                // old object to reallocate
    bool *ok                // true if successful, false otherwise
) ;

void GB_free_memory
(
    void *p,                // pointer to allocated block of memory to free
    size_t nitems,          // number of items to free
    size_t size_of_item     // sizeof each item
) ;

//------------------------------------------------------------------------------
// macros to create/free matrices, vectors, and generic memory
//------------------------------------------------------------------------------

// if PRINT_MALLOC is #define'd above, these macros print diagnostic
// information, meant for development of SuiteSparse:GraphBLAS only

#ifdef PRINT_MALLOC

#define GB_NEW(A,type,vlen,vdim,Ap_option,is_csc,hopt,h,plen)                 \
{                                                                             \
    printf ("\nmatrix new:                   %s = new (%s, vlen = "GBd        \
        ", vdim = "GBd", Ap:%d, csc:%d, hyper:%d %g, plen:"GBd")"             \
        " line %d file %s\n", GB_STR(A), GB_STR(type),                        \
        (int64_t) vlen, (int64_t) vdim, Ap_option, is_csc, hopt, h,           \
        (int64_t) plen, __LINE__, __FILE__) ;                                 \
    info = GB_new (A, type, vlen, vdim, Ap_option, is_csc, hopt, h, plen) ;   \
}

#define GB_CREATE(A,type,vlen,vdim,Ap_option,is_csc,hopt,h,plen,anz,numeric)  \
{                                                                             \
    printf ("\nmatrix create:                %s = new (%s, vlen = "GBd        \
        ", vdim = "GBd", Ap:%d, csc:%d, hyper:%d %g, plen:"GBd", anz:"GBd     \
        " numeric:%d) line %d file %s\n", GB_STR(A), GB_STR(type),            \
        (int64_t) vlen, (int64_t) vdim, Ap_option, is_csc, hopt,              \
        (int64_t) plen, (int64_t) anz, numeric, __LINE__, __FILE__) ;         \
    info = GB_create (A, type, vlen, vdim, Ap_option, is_csc, hopt, h, plen,  \
        anz, numeric) ;                                                       \
}

#define GB_MATRIX_FREE(A)                                                     \
{                                                                             \
    if (A != NULL && *(A) != NULL)                                            \
        printf ("\nmatrix free:                  "                            \
        "matrix_free (%s) line %d file %s\n", GB_STR(A), __LINE__, __FILE__) ;\
    GB_free (A) ;                                                             \
}

#define GB_VECTOR_FREE(v)                                                     \
{                                                                             \
    if (v != NULL && *(v) != NULL)                                            \
        printf ("\nvector free:                  "                            \
        "vector_free (%s) line %d file %s\n", GB_STR(v), __LINE__, __FILE__) ;\
    GB_free ((GrB_Matrix *) v) ;                                              \
}

#define GB_CALLOC_MEMORY(p,n,s)                                               \
    printf ("\ncalloc:                       "                                \
    "%s = calloc (%s = "GBd", %s = "GBd") line %d file %s\n",                 \
    GB_STR(p), GB_STR(n), (int64_t) n, GB_STR(s), (int64_t) s,                \
    __LINE__,__FILE__) ;                                                      \
    p = GB_calloc_memory (n, s) ;

#define GB_MALLOC_MEMORY(p,n,s)                                               \
    printf ("\nmalloc:                       "                                \
    "%s = malloc (%s = "GBd", %s = "GBd") line %d file %s\n",                 \
    GB_STR(p), GB_STR(n), (int64_t) n, GB_STR(s), (int64_t) s,                \
    __LINE__,__FILE__) ;                                                      \
    p = GB_malloc_memory (n, s) ;

#define GB_REALLOC_MEMORY(p,nnew,nold,s,ok)                                   \
{                                                                             \
    printf ("\nrealloc: %14p       "                                          \
    "%s = realloc (%s = "GBd", %s = "GBd", %s = "GBd") line %d file %s\n",    \
    p, GB_STR(p), GB_STR(nnew), (int64_t) nnew, GB_STR(nold), (int64_t) nold, \
    GB_STR(s), (int64_t) s, __LINE__,__FILE__) ;                              \
    p = GB_realloc_memory (nnew, nold, s, p, ok) ;                            \
}

#define GB_FREE_MEMORY(p,n,s)                                                 \
{                                                                             \
    if (p)                                                                    \
    printf ("\nfree:               "                                          \
    "(%s, %s = "GBd", %s = "GBd") line %d file %s\n",                         \
    GB_STR(p), GB_STR(n), (int64_t) n, GB_STR(s), (int64_t) s,                \
    __LINE__,__FILE__) ;                                                      \
    GB_free_memory (p, n, s) ;                                                \
    (p) = NULL ;                                                              \
}

#else

#define GB_NEW(A,type,vlen,vdim,Ap_option,is_csc,hopt,h,plen)                 \
{                                                                             \
    info = GB_new (A, type, vlen, vdim, Ap_option, is_csc, hopt, h, plen) ;   \
}

#define GB_CREATE(A,type,vlen,vdim,Ap_option,is_csc,hopt,h,plen,anz,numeric)  \
{                                                                             \
    info = GB_create (A, type, vlen, vdim, Ap_option, is_csc, hopt, h, plen,  \
        anz, numeric) ;                                                       \
}

#define GB_MATRIX_FREE(A)                                                     \
{                                                                             \
    GB_free (A) ;                                                             \
}

#define GB_VECTOR_FREE(v)                                                     \
{                                                                             \
    GB_free ((GrB_Matrix *) v) ;                                              \
}

#define GB_CALLOC_MEMORY(p,n,s)                                               \
    p = GB_calloc_memory (n, s) ;

#define GB_MALLOC_MEMORY(p,n,s)                                               \
    p = GB_malloc_memory (n, s) ;

#define GB_REALLOC_MEMORY(p,nnew,nold,s,ok)                                   \
{                                                                             \
    p = GB_realloc_memory (nnew, nold, s, p, ok) ;                            \
}

#define GB_FREE_MEMORY(p,n,s)                                                 \
{                                                                             \
    GB_free_memory (p, n, s) ;                                                \
    (p) = NULL ;                                                              \
}

#endif

//------------------------------------------------------------------------------

GrB_Info GB_Entry_print     // print a single value
(
    const GrB_Type type,    // type of value to print
    const void *x           // value to print
) ;

void GB_code_print              // print an entry using a type code
(
    const GB_Type_code code,    // type code of value to print
    const void *x               // entry to print
) ;

GrB_Type GB_code_type           // return the GrB_Type corresponding to the code
(
    const GB_Type_code code,    // type code to convert
    const GrB_Type type         // user type if code is GB_UDT_code
) ;

GrB_Info GB_AxB_alloc           // estimate NNZ(C) and allocate C for C=A*B
(
    GrB_Matrix *Chandle,        // output matrix
    const GrB_Type ctype,       // type of C
    const GrB_Index cvlen,      // vector length of C
    const GrB_Index cvdim,      // # of vectors of C
    const GrB_Matrix M,         // optional mask
    const GrB_Matrix A,         // input matrix A (transposed for dot product)
    const GrB_Matrix B,         // input matrix B
    const bool numeric,         // if true, allocate A->x, else A->x is NULL
    const int64_t rough_guess   // rough estimate of NNZ(C)
) ;

GrB_Info GB_AxB_Gustavson           // C=A*B or C<M>=A*B, Gustavson's method
(
    GrB_Matrix *Chandle,            // output matrix
    const GrB_Matrix M,             // optional matrix
    const GrB_Matrix A_in,          // input matrix A
    const GrB_Matrix B_in,          // input matrix B
    const GrB_Semiring semiring,    // semiring that defines C=A*B
    const bool flipxy               // if true, do z=fmult(b,a) vs fmult(a,b)
//  const bool flops_are_low        // true if flop count is very low
) ;

// used in GB_AxB_heap for temporary workspace
typedef struct
{
    int64_t start ;                 // first entry of A(:,k) is at Ai [start]
    int64_t end ;                   // last entry of A(:,k) is at Ai [end-1]
}
GB_pointer_pair ;

GrB_Info GB_AxB_heap                // C<M>=A*B or C=A*B using a heap
(
    GrB_Matrix *Chandle,            // output matrix
    const GrB_Matrix M,             // mask matrix for C<M>=A*B
    const GrB_Matrix A,             // input matrix
    const GrB_Matrix B,             // input matrix
    const GrB_Semiring semiring,    // semiring that defines C=A*B
    const bool flipxy,              // if true, do z=fmult(b,a) vs fmult(a,b)
    const int64_t bjnz_max          // max # entries in any vector of B
) ;

GrB_Info GB_AxB_saxpy               // C=A*B or C<M>=A*B, gather/scatter or heap
(
    GrB_Matrix *Chandle,            // output matrix
    const GrB_Matrix M,             // optional matrix
    const GrB_Matrix A,             // input matrix A
    const GrB_Matrix B,             // input matrix B
    const GrB_Semiring semiring,    // semiring that defines C=A*B
    const bool flipxy,              // if true, do z=fmult(b,a) vs fmult(a,b)
    const GrB_Desc_Value AxB_method // for auto vs user selection of methods
) ;

GrB_Info GB_AxB_meta                // C<M>=A*B meta algorithm
(
    GrB_Matrix *Chandle,            // output matrix C
    const bool C_is_csc,            // desired CSR/CSC format of C
    GrB_Matrix *MT_handle,          // return MT = M' to caller, if computed
    const GrB_Matrix M_in,          // mask for C<M> (not complemented)
    const GrB_Matrix A_in,          // input matrix
    const GrB_Matrix B_in,          // input matrix
    const GrB_Semiring semiring,    // semiring that defines C=A*B
    bool A_transpose,               // if true, use A', else A
    bool B_transpose,               // if true, use B', else B
    bool flipxy,                    // if true, do z=fmult(b,a) vs fmult(a,b)
    bool *mask_applied,             // if true, mask was applied
    const GrB_Desc_Value AxB_method // for auto vs user selection of methods
) ;

bool GB_semiring_builtin            // true if semiring is builtin
(
    // inputs:
    const GrB_Matrix A,
    const GrB_Matrix B,
    const GrB_Semiring semiring,    // semiring that defines C=A*B
    const bool flipxy,              // true if z=fmult(y,x), flipping x and y
    // outputs, unused by caller if this function returns false
    GB_Opcode *mult_opcode,         // multiply opcode
    GB_Opcode *add_opcode,          // add opcode
    GB_Type_code *xycode,           // type code for x and y inputs
    GB_Type_code *zcode             // type code for z output
) ;

GrB_Info GB_AxB_Gustavson_builtin
(
    GrB_Matrix C,                   // output matrix
    const GrB_Matrix M,             // M matrix for C<M> (not complemented)
    const GrB_Matrix A,             // input matrix
    const GrB_Matrix B,             // input matrix
    const GrB_Semiring semiring,    // semiring that defines C=A*B
    const bool flipxy,              // if true, do z=fmult(b,a) vs fmult(a,b)
    bool *done                      // true if C=A*B has been computed
) ;

GrB_Info GB_AxB_dot                 // C = A'*B using dot product method
(
    GrB_Matrix *Chandle,            // output matrix
    const GrB_Matrix M,             // mask matrix for C<M>=A'*B
    const GrB_Matrix A,             // input matrix
    const GrB_Matrix B,             // input matrix
    const GrB_Semiring semiring,    // semiring that defines C=A*B
    const bool flipxy               // if true, do z=fmult(b,a) vs fmult(a,b)
) ;

GrB_Info GB_mxm                     // C<Mask> = A*B
(
    GrB_Matrix C,                   // input/output matrix for results
    const bool C_replace,           // if true, clear C before writing to it
    const GrB_Matrix Mask,          // optional mask for C, unused if NULL
    const bool Mask_comp,           // if true, use ~Mask
    const GrB_BinaryOp accum,       // optional accum for Z=accum(C,T)
    const GrB_Semiring semiring,    // defines '+' and '*' for C=A*B
    const GrB_Matrix A,             // input matrix
    const bool A_transpose,         // if true, use A' instead of A
    const GrB_Matrix B,             // input matrix
    const bool B_transpose,         // if true, use B' instead of B
    const bool flipxy,              // if true, do z=fmult(b,a) vs fmult(a,b)
    const GrB_Desc_Value AxB_method // for auto vs user selection of methods
) ;

int64_t GB_cumsum               // compute the cumulative sum of an array
(
    int64_t *p,                 // size n+1, undefined on input
    int64_t *count,             // size n+1, input/output
    const int64_t n
) ;

GrB_Info GB_mask                // C<Mask> = Z
(
    GrB_Matrix C,               // both output C and result matrix Cresult
    const GrB_Matrix Mask,      // optional Mask matrix, can be NULL
    GrB_Matrix *Zhandle,        // Z = results of computation, might be shallow
                                // or can even be NULL if Mask is empty and
                                // complemented.  Z is freed when done.
    const bool C_replace,       // true if clear(C) to be done first
    const bool Mask_complement  // true if Mask is to be complemented
) ;

GrB_Info GB_accum_mask          // C<M> = accum (C,T)
(
    GrB_Matrix C,               // input/output matrix for results
    const GrB_Matrix M_in,      // optional mask for C, unused if NULL
    const GrB_Matrix MT_in,     // MT=M' if computed already in the caller
    const GrB_BinaryOp accum,   // optional accum for Z=accum(C,results)
    GrB_Matrix *Thandle,        // results of computation, freed when done
    const bool C_replace,       // if true, clear C first
    const bool Mask_complement  // if true, complement the mask
) ;

GrB_Info GB_Descriptor_get      // get the contents of a descriptor
(
    const GrB_Descriptor desc,  // descriptor to query, may be NULL
    bool *C_replace,            // if true replace C before C<Mask>=Z
    bool *Mask_comp,            // if true use logical negation of Mask
    bool *In0_transpose,        // if true transpose first input
    bool *In1_transpose,        // if true transpose second input
    GrB_Desc_Value *AxB_method  // method for C=A*B
) ;

GrB_Info GB_compatible          // SUCCESS if all is OK, *_MISMATCH otherwise
(
    const GrB_Type ctype,       // the type of C (matrix or scalar)
    const GrB_Matrix C,         // the output matrix C; NULL if C is a scalar
    const GrB_Matrix Mask,      // optional Mask, NULL if no mask
    const GrB_BinaryOp accum,   // C<Mask> = accum(C,T) is computed
    const GrB_Type ttype        // type of T
) ;

GrB_Info GB_Mask_compatible     // check type and dimensions of mask
(
    const GrB_Matrix Mask,      // mask to check
    const GrB_Matrix C,         // C<Mask>= ...
    const GrB_Index nrows,      // size of output if C is NULL
    const GrB_Index ncols
) ;

GrB_Info GB_BinaryOp_compatible     // check for domain mismatch
(
    const GrB_BinaryOp op,          // binary operator to check
    const GrB_Type ctype,           // C must be compatible with op->ztype
    const GrB_Type atype,           // A must be compatible with op->xtype
    const GrB_Type btype,           // B must be compatible with op->ytype
    const GB_Type_code bcode        // B may not have a type, just a code
) ;

void GB_qsort_1         // sort array A of size 1-by-n
(
    int64_t A_0 [ ],    // size-n array
    const int64_t n
) ;

void GB_qsort_2a        // sort array A of size 2-by-n, using 1 key (A [0][])
(
    int64_t A_0 [ ],    // size n array
    int64_t A_1 [ ],    // size n array
    const int64_t n
) ;

void GB_qsort_2b        // sort array A of size 2-by-n, using 2 keys (A [0:1][])
(
    int64_t A_0 [ ],    // size n array
    int64_t A_1 [ ],    // size n array
    const int64_t n
) ;

void GB_qsort_3         // sort array A of size 3-by-n, using 3 keys (A [0:2][])
(
    int64_t A_0 [ ],    // size n array
    int64_t A_1 [ ],    // size n array
    int64_t A_2 [ ],    // size n array
    const int64_t n
) ;

GrB_Info GB_subref_numeric      // C = A (I,J), extract the values
(
    GrB_Matrix *Chandle,        // output C
    const bool C_is_csc,        // requested format of C
    const GrB_Matrix A,         // input matrix
    const GrB_Index *I,         // list of indices, duplicates OK
    int64_t ni,                 // length of I
    const GrB_Index *J,         // list of vector indices, duplicates OK
    int64_t nj                  // length of J
    , const bool must_sort      // if true, must return C sorted
) ;

GrB_Info GB_subref_symbolic     // C = A (I,J), extract the pattern
(
    GrB_Matrix *Chandle,        // output C
    const bool C_is_csc,        // requested format of C
    const GrB_Matrix A,         // input matrix
    const GrB_Index *I,         // list of indices, duplicates OK
    int64_t ni,                 // length of I
    const GrB_Index *J,         // list of vector indices, duplicates OK
    int64_t nj                  // length of J
) ;

GrB_Info GB_extract                 // C<Mask> = accum (C, A(I,J))
(
    GrB_Matrix C,                   // input/output matrix for results
    const bool C_replace,           // C matrix descriptor
    const GrB_Matrix Mask,          // optional mask for C, unused if NULL
    const bool Mask_comp,           // Mask descriptor
    const GrB_BinaryOp accum,       // optional accum for Z=accum(C,T)
    const GrB_Matrix A,             // input matrix
    const bool A_transpose,         // A matrix descriptor
    const GrB_Index *Rows,          // row indices
    const GrB_Index nRows_in,       // number of row indices
    const GrB_Index *Cols,          // column indices
    const GrB_Index nCols_in        // number of column indices
) ;

GrB_Info GB_eWise                   // C<M> = accum (C, A+B) or A.*B
(
    GrB_Matrix C,                   // input/output matrix for results
    const bool C_replace,           // if true, clear C before writing to it
    const GrB_Matrix M,             // optional mask for C, unused if NULL
    const bool Mask_comp,           // if true, complement the mask M
    const GrB_BinaryOp accum,       // optional accum for Z=accum(C,T)
    const GrB_BinaryOp op,          // defines '+' for C=A+B, or .* for A.*B
    const GrB_Matrix A,             // input matrix
    bool A_transpose,               // if true, use A' instead of A
    const GrB_Matrix B,             // input matrix
    bool B_transpose,               // if true, use B' instead of B
    const bool eWiseAdd             // if true, do set union (like A+B),
                                    // otherwise do intersection (like A.*B)
) ;

GrB_Info GB_reduce_to_column        // w<mask> = accum (w,reduce(A))
(
    GrB_Matrix w,                   // input/output for results, size n-by-1
    const GrB_Matrix mask,          // optional mask for w, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for z=accum(w,t)
    const GrB_BinaryOp reduce,      // reduce operator for t=reduce(A)
    const GrB_Matrix A,             // first input:  matrix A
    const GrB_Descriptor desc       // descriptor for w, mask, and A
) ;

GrB_Info GB_reduce_to_scalar    // t = reduce_to_scalar (A)
(
    void *c,                    // result scalar
    const GrB_Type ctype,       // the type of the scalar c
    const GrB_BinaryOp accum,   // optional accum for c = accum(c,t)
    const GrB_Monoid reduce,    // monoid to do the reduction
    const GrB_Matrix A          // matrix to reduce
) ;

GB_Opcode GB_boolean_rename     // renamed opcode
(
    const GB_Opcode opcode      // opcode to rename
) ;

bool GB_Index_multiply      // true if ok, false if overflow
(
    GrB_Index *c,           // c = a*b, or zero if overflow occurs
    const int64_t a,
    const int64_t b
) ;

bool GB_size_t_multiply     // true if ok, false if overflow
(
    size_t *c,              // c = a*b, or zero if overflow occurs
    const size_t a,
    const size_t b
) ;

GrB_Info GB_extractTuples       // extract all tuples from a matrix
(
    GrB_Index *I,               // array for returning row indices of tuples
    GrB_Index *J,               // array for returning col indices of tuples
    void *X,                    // array for returning values of tuples
    GrB_Index *nvals,           // I,J,X size on input; # tuples on output
    const GB_Type_code xcode,   // type of array X
    const GrB_Matrix A          // matrix to extract tuples from
) ;

GrB_Info GB_extractElement      // extract a single entry, x = A(i,j)
(
    void *x,                    // scalar to extract, not modified if not found
    const GB_Type_code xcode,   // type of the scalar x
    const GrB_Matrix A,         // matrix to extract a scalar from
    const GrB_Index i,          // row index
    const GrB_Index j           // column index
) ;

GrB_Info GB_Monoid_new          // create a monoid
(
    GrB_Monoid *monoid_handle,  // handle of monoid to create
    const GrB_BinaryOp op,      // binary operator of the monoid
    const void *id,             // identity value
    const GB_Type_code idcode   // identity code
) ;

GrB_Info GB_user_build          // check inputs then build matrix
(
    GrB_Matrix C,               // matrix to build
    const GrB_Index *I,         // row indices of tuples
    const GrB_Index *J,         // col indices of tuples
    const void *S,              // array of values of tuples
    const GrB_Index nvals,      // number of tuples
    const GrB_BinaryOp dup,     // binary function to assemble duplicates
    const GB_Type_code scode,   // GB_Type_code of S array
    const bool is_matrix        // true if C is a matrix, false if GrB_Vector
) ;

GrB_Info GB_build               // build matrix
(
    GrB_Matrix C,               // matrix to build
    const GrB_Index *I_in,      // row indices of tuples
    const GrB_Index *J_in,      // col indices of tuples
    const void *X,              // array of values of tuples
    const GrB_Index nvals,      // number of tuples
    const GrB_BinaryOp dup,     // binary function to assemble duplicates
    const GB_Type_code X_code,  // GB_Type_code of X array
    const bool is_matrix,       // true if C is a matrix, false if GrB_Vector
    const bool ijcheck          // true if I and J are to be checked
) ;

GrB_Info GB_builder
(
    GrB_Matrix *Thandle,            // matrix T to build
    const GrB_Type ttype,           // type of output matrix T
    const int64_t vlen,             // length of each vector of T
    const int64_t vdim,             // number of vectors in T
    const bool is_csc,              // true if T is CSC, false if CSR
    int64_t **iwork_handle,         // for (i,k) or (j,i,k) tuples
    int64_t **jwork_handle,         // for (j,i,k) tuples
    const bool already_sorted,      // true if tuples already sorted on input
    const void *S,                  // array of values of tuples
    const int64_t len,              // number of tuples, and size of kwork
    const int64_t ijlen,            // size of iwork and jwork arrays
    const GrB_BinaryOp dup,         // binary function to assemble duplicates,
                                    // if NULL use the "SECOND" function to
                                    // keep the most recent duplicate.
    const GB_Type_code scode        // GB_Type_code of S array
) ;

GrB_Info GB_build_factory           // build a matrix
(
    GrB_Matrix *Thandle,            // matrix T to build
    const int64_t tnz0,             // final NNZ(T)
    int64_t **iwork_handle,         // for (i,k) or (j,i,k) tuples
    int64_t **kwork_handle,         // for (i,k) or (j,i,k) tuples
    const void *S,                  // array of values of tuples
    const int64_t len,              // number of tuples and size of kwork
    const int64_t ijlen,            // size of iwork array
    const GrB_BinaryOp dup,         // binary function to assemble duplicates,
                                    // if NULL use the "SECOND" function to
                                    // keep the most recent duplicate.
    const GB_Type_code scode        // GB_Type_code of S array
) ;

GrB_Info GB_wait                // finish all pending computations
(
    GrB_Matrix A                // matrix with pending computations
) ;

GrB_Info GB_pending_add             // add a pending tuple A(i,j) to a matrix
(
    GrB_Matrix A,                   // matrix to modify
    const void *scalar,             // scalar to add to the pending tuples of A
    const GrB_Type stype,           // scalar type
    const GrB_BinaryOp pop,         // new A->operator_pending, if 1st pending
    const int64_t i,                // index into vector
    const int64_t j                 // vector index
) ;

void GB_pending_free            // free all pending tuples
(
    GrB_Matrix A                // matrix with pending tuples to free
) ;

void GB_queue_remove            // remove matrix from queue
(
    GrB_Matrix A                // matrix to remove
) ;

void GB_queue_insert            // insert matrix at the head of queue
(
    GrB_Matrix A                // matrix to insert
) ;

GrB_Matrix GB_queue_remove_head ( ) ; // return matrix or NULL if queue empty

void GB_queue_check
(
    GrB_Matrix A,           // matrix to check
    GrB_Matrix *head,       // head of the queue
    GrB_Matrix *prev,       // prev from A
    GrB_Matrix *next,       // next after A
    bool *enqd              // true if A is in the queue
) ;

GrB_Info GB_setElement          // set a single entry, C(i,j) = x
(
    GrB_Matrix C,               // matrix to modify
    const void *x,              // scalar to set
    const GrB_Index i_in,       // row index
    const GrB_Index j_in,       // column index
    const GB_Type_code xcode    // type of the scalar x
) ;

GrB_Info GB_block   // apply all pending computations if blocking mode enabled
(
    GrB_Matrix A
) ;

GrB_Info GB_subassign               // C(Rows,Cols)<M> += A or A'
(
    GrB_Matrix C,                   // input/output matrix for results
    const bool C_replace,           // descriptor for C
    const GrB_Matrix M_in,          // optional mask for C(Rows,Cols)
    const bool Mask_comp,           // true if mask is complemented
    bool M_transpose,               // true if the mask should be transposed
    const GrB_BinaryOp accum,       // optional accum for accum(C,T)
    const GrB_Matrix A_in,          // input matrix
    bool A_transpose,               // true if A is transposed
    const GrB_Index *Rows,          // row indices
    const GrB_Index nRows_in,       // number of row indices
    const GrB_Index *Cols,          // column indices
    const GrB_Index nCols_in,       // number of column indices
    const bool scalar_expansion,    // if true, expand scalar to A
    const void *scalar,             // scalar to be expanded
    const GB_Type_code scalar_code  // type code of scalar to expand
) ;

GrB_Info GB_assign                  // C<M>(Rows,Cols) += A or A'
(
    GrB_Matrix C,                   // input/output matrix for results
    const bool C_replace,           // descriptor for C
    const GrB_Matrix M_in,          // optional mask for C
    const bool Mask_comp,           // true if mask is complemented
    bool M_transpose,               // true if the mask should be transposed
    const GrB_BinaryOp accum,       // optional accum for accum(C,T)
    const GrB_Matrix A_in,          // input matrix
    bool A_transpose,               // true if A is transposed
    const GrB_Index *Rows,          // row indices
    const GrB_Index nRows_in,       // number of row indices
    const GrB_Index *Cols,          // column indices
    const GrB_Index nCols_in,       // number of column indices
    const bool scalar_expansion,    // if true, expand scalar to A
    const void *scalar,             // scalar to be expanded
    const GB_Type_code scalar_code  // type code of scalar to expand
    , const bool col_assign         // true for GrB_Col_assign
    , const bool row_assign         // true for GrB_Row_assign
) ;

GrB_Info GB_subassign_kernel        // C(I,J)<M> = A or accum (C (I,J), A)
(
    GrB_Matrix C,                   // input/output matrix for results
    bool C_replace,                 // C matrix descriptor
    const GrB_Matrix M,             // optional mask for C(I,J), unused if NULL
    const bool Mask_comp,           // mask descriptor
    const GrB_BinaryOp accum,       // optional accum for Z=accum(C(I,J),A)
    const GrB_Matrix A,             // input matrix (NULL for scalar expansion)
    const GrB_Index *I,             // list of indices
    const int64_t   ni,             // number of indices
    const GrB_Index *J,             // list of vector indices
    const int64_t   nj,             // number of column indices
    const bool scalar_expansion,    // if true, expand scalar to A
    const void *scalar,             // scalar to be expanded
    const GB_Type_code scalar_code  // type code of scalar to expand
) ;

GrB_Info GB_subassign_scalar        // C(I,J)<Mask> = accum (C(I,J),x)
(
    GrB_Matrix C,                   // input/output matrix for results
    const GrB_Matrix Mask,          // optional mask for C(I,J), unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(C(I,J),T)
    const void *scalar,             // scalar to assign to C(I,J)
    const GB_Type_code scalar_code, // type code of scalar to assign
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Index *J,             // column indices
    const GrB_Index nj,             // number of column indices
    const GrB_Descriptor desc       // descriptor for C(I,J) and Mask
) ;

GrB_Info GB_assign_scalar           // C<Mask>(I,J) = accum (C(I,J),x)
(
    GrB_Matrix C,                   // input/output matrix for results
    const GrB_Matrix Mask,          // optional mask for C, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(C(I,J),T)
    const void *scalar,             // scalar to assign to C(I,J)
    const GB_Type_code scalar_code, // type code of scalar to assign
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Index *J,             // column indices
    const GrB_Index nj,             // number of column indices
    const GrB_Descriptor desc       // descriptor for C and Mask
) ;

bool GB_op_is_second    // return true if op is SECOND, of the right type
(
    GrB_BinaryOp op,
    GrB_Type type
) ;


//------------------------------------------------------------------------------

char *GB_code_string            // return a static string for a type name
(
    const GB_Type_code code     // code to convert to string
) ;

GrB_Info GB_resize              // change the size of a matrix
(
    GrB_Matrix A,               // matrix to modify
    const GrB_Index nrows_new,  // new number of rows in matrix
    const GrB_Index ncols_new   // new number of columns in matrix
) ;

GrB_Info GB_kron                    // C<Mask> = accum (C, kron(A,B))
(
    GrB_Matrix C,                   // input/output matrix for results
    const bool C_replace,           // if true, clear C before writing to it
    const GrB_Matrix Mask,          // optional mask for C, unused if NULL
    const bool Mask_comp,           // if true, use ~Mask
    const GrB_BinaryOp accum,       // optional accum for Z=accum(C,T)
    const GrB_BinaryOp op,          // defines '*' for kron(A,B)
    const GrB_Matrix A,             // input matrix
    bool A_transpose,               // if true, use A' instead of A
    const GrB_Matrix B,             // input matrix
    bool B_transpose                // if true, use B' instead of B
) ;

GrB_Info GB_kron_kernel             // C = kron (A,B)
(
    GrB_Matrix *Chandle,            // output matrix
    const bool C_is_csc,            // desired format of C
    const GrB_BinaryOp op,          // multiply operator
    const GrB_Matrix A,             // input matrix
    const GrB_Matrix B              // input matrix
) ;

void GB_stats
(
    GxB_Statistics *stats           // return statistics
) ;

void GB_apply_op            // apply a unary operator, Cx = op ((xtype) Ax)
(
    void *Cx,               // output array, of type op->ztype
    const GrB_UnaryOp op,   // operator to apply
    const void *Ax,         // input array, of type atype
    const GrB_Type atype,   // type of Ax
    const int64_t anz       // size of Ax and Cx
) ;

GrB_Info GB_transpose           // C=A', C=(ctype)A or C=op(A')
(
    GrB_Matrix *Chandle,        // output matrix C, possibly modified in place
    GrB_Type ctype,             // desired type of C; if NULL use A->type.
                                // ignored if op is present (cast to op->ztype)
    const bool C_is_csc,        // desired CSR/CSC format of C
    const GrB_Matrix A_in,      // input matrix, NULL if C is modified in place
    const GrB_UnaryOp op_in     // optional operator to apply to the values
) ;

int64_t GB_nvec_nonempty        // return # of non-empty vectors
(
    const GrB_Matrix A          // input matrix to examine
) ;

GrB_Info GB_to_nonhyper     // convert a matrix to non-hypersparse
(
    GrB_Matrix A            // matrix to convert to non-hypersparse
) ;

GrB_Info GB_to_hyper        // convert a matrix to hypersparse
(
    GrB_Matrix A            // matrix to convert to hypersparse
) ;

bool GB_to_nonhyper_check   // check for conversion to hypersparse
(
    GrB_Matrix A,           // matrix to check
    int64_t k,              // # of non-empty vectors of A, an estimate is OK,
                            // but normally A->nvec_nonempty
    int64_t vdim            // normally A->vdim
) ;

bool GB_to_hyper_check      // check for conversion to hypersparse
(
    GrB_Matrix A,           // matrix to check
    int64_t k,              // # of non-empty vectors of A, an estimate is OK,
                            // but normally A->nvec_nonempty
    int64_t vdim            // normally A->vdim
) ;

GrB_Info GB_to_hyper_conform    // conform a matrix to its desired format
(
    GrB_Matrix A                // matrix to conform
) ;

void GB_global_option_set       // set one or more global options
(
    bool set_hyper_ratio,       // if true, set the hyper_ratio
    double hyper_ratio,
    bool set_format,            // if true, set the format
    GxB_Format_Value format
) ;

void GB_global_option_get       // get one or more global options
(
    double *hyper_ratio,
    GxB_Format_Value *format,
    bool *is_csc
) ;

//------------------------------------------------------------------------------
// Global storage: for all threads in a user application that uses GraphBLAS
//------------------------------------------------------------------------------

// Global storage is used to record a list of matrices with pending operations
// (fo GrB_wait), and to keep track of the GraphBLAS mode (blocking or
// non-blocking).

typedef struct
{

    //--------------------------------------------------------------------------
    // queue of matrices with work to do
    //--------------------------------------------------------------------------

    // In non-blocking mode, GraphBLAS needs to keep track of all matrices that
    // have pending operations that have not yet been finished.  In the current
    // implementation, these are matrices with pending tuples from
    // GrB_setElement, GxB_subassign, and GrB_assign that haven't been added to
    // the matrix yet.

    // A matrix with no pending tuples is not in the list.  When a matrix gets
    // its first pending tuple, it is added to the list.  A matrix is removed
    // from the list when another operation needs to use the matrix; in that
    // case the pending tuples are assembled for just that one matrix.  The
    // GrB_wait operation iterates through the entire list and assembles all
    // the pending tuples for all the matrices in the list, leaving the list
    // emtpy.  A simple link list suffices for the list.  The links are in the
    // matrices themselves so no additional memory needs to be allocated.  The
    // list never needs to be searched; if a particular matrix is to be removed
    // from the list, GraphBLAS has already been given the matrix handle, and
    // the prev & next pointers it contains.  All of these operations can thus
    // be done in O(1) time, except for GrB_wait which needs to traverse the
    // whole list once and then the list is empty afterwards.

    void *queue_head ;          // head pointer to matrix queue

    GrB_Mode mode ;             // GrB_NONBLOCKING or GrB_BLOCKING

    bool GrB_init_called ;      // true if GrB_init already called

    //--------------------------------------------------------------------------
    // hypersparsity control
    //--------------------------------------------------------------------------

    // See GB_matrix.h and GraphBLAS.h for a description.

    double hyper_ratio ;        // default hyper_ratio for new matrices
    bool is_csc ;               // default CSR/CSC format for new matrices

    //--------------------------------------------------------------------------
    // malloc tracking
    //--------------------------------------------------------------------------

    // nmalloc:  To aid in searching for memory leaks, GraphBLAS keeps track of
    // the number of blocks of allocated by malloc, calloc, or realloc that
    // have not yet been freed.  The count starts at zero.  malloc and calloc
    // increment this count, and free (of a non-NULL pointer) decrements it.
    // realloc increments the count it if is allocating a new block, but it
    // does this by calling GB_malloc_memory.

    // inuse: the # of bytes currently in use by all threads

    // maxused: the max value of inuse since the call to GrB_init

    // malloc_debug: this is used for testing only (GraphBLAS/Tcov).  If true,
    // then use malloc_debug_count for testing memory allocation and
    // out-of-memory conditions.  If malloc_debug_count > 0, the value is
    // decremented after each allocation of memory.  If malloc_debug_count <=
    // 0, the GB_*_memory routines pretend to fail; returning NULL and not
    // allocating anything.

    int64_t nmalloc ;               // number of blocks allocated but not freed
    bool malloc_debug ;             // if true, test memory hanlding
    int64_t malloc_debug_count ;    // for testing memory handling
    int64_t inuse ;                 // memory space current in use
    int64_t maxused ;               // high water memory usage

}
GB_Global_struct ;

extern GB_Global_struct GB_Global ;

//------------------------------------------------------------------------------
// Thread local storage
//------------------------------------------------------------------------------

// Thread local storage is used to to record the details of the last error
// encountered (for GrB_error).  If the user application is multi-threaded,
// each thread that calls GraphBLAS needs its own private copy of these
// variables.

#define GB_RLEN 3000
#define GB_DLEN 2048

typedef struct
{

    //--------------------------------------------------------------------------
    // Error status
    //--------------------------------------------------------------------------

    GrB_Info info ;             // last error code
    GrB_Index row, col ;        // last A(row,col) searched for but not found
    bool is_matrix ;            // last search matrix (true) or vector (false)
    const char *where ;         // GraphBLAS function where error occurred
    const char *file ;          // GraphBLAS filename where error occured
    int line ;                  // line in the GraphBLAS file of error
    char details [GB_DLEN+1] ;  // string with details of the error
    char report [GB_RLEN+1] ;   // string returned by GrB_error

    //--------------------------------------------------------------------------
    // workspace
    //--------------------------------------------------------------------------

    // Initialized space: Mark is an array of size Mark_size.  When cleared,
    // Mark [0..Mark_size-1] < Mark_flag holds.

    int64_t *Mark ;                 // initialized space
    int64_t Mark_flag ;             // current watermark in Mark [...]
    int64_t Mark_size ;             // current size of Mark array

    // Uninitialized space:
    void *Work ;                    // uninitialized space
    size_t Work_size ;              // current size of Work array

    // Initialized space: Flag is an array of size Flag_size.  When cleared,
    // it is all equal to zero.

    int8_t *Flag ;                  // initialized space
    int64_t Flag_size ;             // current size of Flag array

    //--------------------------------------------------------------------------
    // random seed for GB_rand
    //--------------------------------------------------------------------------

    uint64_t seed ;

    //--------------------------------------------------------------------------
    // algorithm selector
    //--------------------------------------------------------------------------

    GrB_Desc_Value AxB_method ;     // last method used in C=A*B

}
GB_thread_local_struct ;

_Thread_local extern GB_thread_local_struct GB_thread_local ;

//------------------------------------------------------------------------------
// random number generator
//------------------------------------------------------------------------------

// return a random GrB_Index, in range 0 to 2^60
#define GB_RAND_MAX 32767

// return a random number between 0 and GB_RAND_MAX
static inline GrB_Index GB_rand15 ( )
{ 
   GB_thread_local.seed = GB_thread_local.seed * 1103515245 + 12345 ;
   return ((GB_thread_local.seed / 65536) % (GB_RAND_MAX + 1)) ;
}

// return a random GrB_Index, in range 0 to 2^60
static inline GrB_Index GB_rand ( )
{ 
    GrB_Index i = GB_rand15 ( ) ;
    i = GB_RAND_MAX * i + GB_rand15 ( ) ;
    i = GB_RAND_MAX * i + GB_rand15 ( ) ;
    i = GB_RAND_MAX * i + GB_rand15 ( ) ;
    return (i) ;
}

//------------------------------------------------------------------------------
// error logging
//------------------------------------------------------------------------------

// The ERROR and LOG macros work together.  If an error occurs, the ERROR macro
// records the details in thread local storage (GB_thread_local), and returns
// the GrB_info to its 'caller'.  This value can then be returned, or set to
// an info variable of type GrB_Info.  For example:
//
//  if (i >= nrows)
//  {
//      return (ERROR (GrB_INDEX_OUT_OF_BOUNDS, (LOG,
//          "Row index %d out of bounds; must be < %d", i, nrows))) ;
//  }
//
// The user can then do:
//
//  printf ("%s", GrB_error ( )) ;
//
// To print details of the error, which includes: which user-callable function
// encountered the error, the error status (GrB_INDEX_OUT_OF_BOUNDS), the
// details ("Row index 102 out of bounds, must be < 100"), and finally the
// exact GraphBLAS filename and line number where the error was caught.

#define LOG GB_thread_local.details, GB_DLEN
#define ERROR(f,s)                          \
(                                           \
        snprintf s ,                        \
        (GB_thread_local.file = __FILE__ ), \
        (GB_thread_local.line = __LINE__ ), \
        (GB_thread_local.info = f)          \
)

// return (OUT_OF_MEMORY (mem)) ; reports an out-of-memory error, including
// the amount of memory needed for this phase of computation
#define OUT_OF_MEMORY(mem) ERROR (GrB_OUT_OF_MEMORY, (LOG,  \
    "out of memory, %g GBytes required", (double) (mem)))

// return (NO_MEMORY) ; reports an out-of-memory error, but does not pring
// the memory needed (normally very small)
#define NO_MEMORY ERROR (GrB_OUT_OF_MEMORY, (LOG, "out of memory"))

// The WHERE macro keeps track of the currently running user-callable function.
// User-callable functions in this implementation are written so that they do
// not call other unrelated user-callable functions (except for GrB_*free).
// Related user-callable functions can call each other since they all report
// the same type-generic name.  Internal functions can be called by many
// different user-callable functions, directly or indirectly.  It would not be
// helpful to report the name of an internal function that flagged an error
// condition.  Thus, each time a user-callable function is entered (except
// GrB_*free), it logs the name of the function with the WHERE macro.
// GrB_*free does not encounter error conditions so it doesn't need to be
// logged by the WHERE macro.

#define WHERE(w) { GB_thread_local.where = w ; }

// The REPORT_NO_VALUE macro is used only by GB_extractElement, to report
// the last A(i,j) that was searched but which was not present in the matrix.
// In the interest of speed, no other action is taken.   If A(i,j) is found,
// then the status is set to GrB_SUCCESS via the REPORT_SUCCESS macro.

#define REPORT_NO_VALUE(i,j)                    \
(                                               \
    GB_thread_local.row = i,                    \
    GB_thread_local.col = j,                    \
    GB_thread_local.info = GrB_NO_VALUE         \
)

// REPORT_MATRIX and REPORT_VECTOR is only used in GrB_*_extractElement_*
// to allow the error reporting to distinguish between the GrB_Matrix_* and
// GrB_Vector_* versions of the functions.
#define REPORT_MATRIX(info)                                             \
    if (info != GrB_SUCCESS) GB_thread_local.is_matrix = true ;         \

#define REPORT_VECTOR(info)                                             \
    if (info != GrB_SUCCESS) GB_thread_local.is_matrix = false ;        \

// All methods and operations use this macro for a successful return
#define REPORT_SUCCESS                          \
(                                               \
    GB_thread_local.info = GrB_SUCCESS          \
)

//------------------------------------------------------------------------------
// Thread local storage management
//------------------------------------------------------------------------------

// All allocator and free functions are called *_walloc and *_wfree.

GrB_Info GB_Mark_walloc             // allocate Mark space
(
    int64_t Mark_required           // ensure Mark is at least this large
) ;

void GB_Mark_wfree ( ) ;            // free the Mark array

int64_t GB_Mark_reset   // increment the Mark_flag by reset
(
    int64_t reset,      // does Mark_flag += reset
    int64_t range       // clear Mark if flag+reset+range overflows
) ;

#ifndef NDEBUG
    #define ASSERT_MARK_IS_RESET                                            \
    {                                                                       \
        for (int64_t i = 0 ; i < GB_thread_local.Mark_size ; i++)           \
        {                                                                   \
            ASSERT (GB_thread_local.Mark [i] < GB_thread_local.Mark_flag) ; \
        }                                                                   \
    }
#else
    #define ASSERT_MARK_IS_RESET
#endif

GrB_Info GB_Work_walloc             // allocate Work space
(
    size_t Work_nitems_required,    // ensure Work is at least this large,
    size_t Work_itemsize            // # items times size of each item
) ;

void GB_Work_wfree ( ) ;            // free the Work array

GrB_Info GB_Flag_walloc             // allocate Flag space
(
    int64_t Flag_required           // ensure Flag is at least this large
) ;

void GB_Flag_wfree ( ) ;            // free the Flag array

#ifndef NDEBUG
    #define ASSERT_FLAG_IS_CLEAR                                            \
    {                                                                       \
        for (int64_t i = 0 ; i < GB_thread_local.Flag_size ; i++)           \
        {                                                                   \
            ASSERT (GB_thread_local.Flag [i] == 0) ;                        \
        }                                                                   \
    }
#else
    #define ASSERT_FLAG_IS_CLEAR
#endif

void GB_wfree ( ) ;                 // free all thread-local storage

//------------------------------------------------------------------------------
// boiler plate macros for checking inputs and returning if an error occurs
//------------------------------------------------------------------------------

// Functions use these macros to check/get their inputs and return an error
// if something has gone wrong.

// check if a required arg is NULL
#define RETURN_IF_NULL(arg)                                             \
    if ((arg) == NULL)                                                  \
    {                                                                   \
        /* the required arg is NULL */                                  \
        return (ERROR (GrB_NULL_POINTER, (LOG,                          \
            "Required argument is null: [%s]", GB_STR(arg)))) ;         \
    }

// arg may be NULL, but if non-NULL then it must be initialized
#define RETURN_IF_FAULTY(arg)                                           \
    if ((arg) != NULL && (arg)->magic != MAGIC)                         \
    {                                                                   \
        if ((arg)->magic == MAGIC2)                                     \
        {                                                               \
            /* optional arg is not NULL, but invalid */                 \
            return (ERROR (GrB_INVALID_OBJECT, (LOG,                    \
                "Argument is invalid: [%s]", GB_STR(arg)))) ;           \
        }                                                               \
        else                                                            \
        {                                                               \
            /* optional arg is not NULL, but not initialized */         \
            return (ERROR (GrB_UNINITIALIZED_OBJECT, (LOG,              \
                "Argument is uninitialized: [%s]", GB_STR(arg)))) ;     \
        }                                                               \
    }

// arg must not be NULL, and it must be initialized
#define RETURN_IF_NULL_OR_FAULTY(arg)                                   \
    RETURN_IF_NULL (arg) ;                                              \
    RETURN_IF_FAULTY (arg) ;                                            \

// check the descriptor and extract its contents
#define GET_DESCRIPTOR(info,desc,dout,dm,d0,d1,dalgo)                   \
    GrB_Info info ;                                                     \
    bool dout, dm, d0, d1 ;                                             \
    GrB_Desc_Value dalgo ;                                              \
    /* if desc is NULL then defaults are used.  This is OK */           \
    info = GB_Descriptor_get (desc, &dout, &dm, &d0, &d1, &dalgo) ;     \
    if (info != GrB_SUCCESS)                                            \
    {                                                                   \
        /* desc not NULL, but uninitialized or an invalid object */     \
        return (info) ;                                                 \
    }

// C<Mask>=Z ignores Z if an empty Mask is complemented, so return from
// the method without computing anything.  But do apply the mask.
#define RETURN_IF_QUICK_MASK(C, C_replace, Mask, Mask_comp)             \
    if (Mask_comp && Mask == NULL)                                      \
    {                                                                   \
        /* C<~NULL>=NULL since result does not depend on computing Z */ \
        return (GB_mask (C, NULL, NULL, C_replace, true)) ;             \
    }

//------------------------------------------------------------------------------
// Pending upddate and zombies
//------------------------------------------------------------------------------

// FLIP is a kind of  "negation" about (-1) of a zero-based index.
// FLIP(i) acts like -i for one-based indices.
// If i >= 0 then it is not flipped.
// If i < 0 then it has been flipped.
// Like negation, FLIP is its own inverse: FLIP (FLIP (i)) == i.
// The "nil" value, -1, doesn't change when flipped: FLIP (-1) = -1.
// UNFLIP(i) is like taking an absolute value, undoing any flip.

// An entry A(i,j) in a matrix can be marked as a "zombie".  A zombie is an
// entry that has been marked for deletion, but hasn't been deleted yet because
// it's more efficient to delete all zombies all at once, instead of one at a
// time.  Zombies are created by submatrix assignment, C(I,J)=A which copies
// not only new entries into C, but it also deletes entries already present in
// C.  If an entry appears in A but not C(I,J), it is a new entry; new entries
// placed in the pending tuple lists to be added later.  If an entry appear in
// C(I,J) but NOT in A, then it is marked for deletion by flipping its row
// index, marking it as a zombie.

// Zombies can be restored as regular entries by GrB_*assign.  If an assignment
// C(I,J)=A finds an entry in A that is a zombie in C, the zombie becomes a
// regular entry, taking on the value from A.  The row index is unflipped.

// Zombies are deleted and pending tuples are added into the matrix all at
// once, by GB_wait.

#define NIL                 (-1)
#define FLIP(i)             (-(i)-2)
#define IS_FLIPPED(i)       ((i) < 0)
#define IS_ZOMBIE(i)        ((i) < 0)
#define IS_NOT_FLIPPED(i)   ((i) >= 0)
#define IS_NOT_ZOMBIE(i)    ((i) >= 0)
#define UNFLIP(i)           (((i) < 0) ? FLIP(i) : (i))

// true if a matrix has pending tuples
#define PENDING(A) ((A) != NULL && (A)->n_pending > 0)

// true if a matrix is allowed to have pending tuples
#define PENDING_OK(A) ((A) != NULL && (A)->n_pending >= 0)

// true if a matrix has zombies
#define ZOMBIES(A) ((A) != NULL && (A)->nzombies > 0)

// true if a matrix is allowed to have zombies
#define ZOMBIES_OK(A) ((A) != NULL && (A)->nzombies >= 0)

// do all pending updates:  delete zombies and assemble any pending tuples
#define GB_WAIT(A)                                                      \
{                                                                       \
    GrB_Info info = GB_wait ((GrB_Matrix) A) ;                          \
    if (info != GrB_SUCCESS)                                            \
    {                                                                   \
        /* out of memory; no other error possible */                    \
        ASSERT (info == GrB_OUT_OF_MEMORY) ;                            \
        return (info) ;                                                 \
    }                                                                   \
    ASSERT (!ZOMBIES (A)) ;                                             \
    ASSERT (!PENDING (A)) ;                                             \
}

// wait for any pending operations: both pending tuples and zombies
#define WAIT(A)                                                         \
{                                                                       \
    if (PENDING (A) || ZOMBIES (A)) GB_WAIT (A) ;                       \
}

// just wait for pending tuples; zombies are OK but removed anyway if the
// matrix also has pending tuples.  They are left in place if the matrix has
// zombies but no pending tuples.
#define WAIT_PENDING(A)                                                 \
{                                                                       \
    if (PENDING (A)) GB_WAIT (A) ;                                      \
    ASSERT (ZOMBIES_OK (A)) ;                                           \
}

// true if a matrix has no entries; zombies OK
#define EMPTY(A) ((NNZ (A) == 0) && !PENDING (A))

//------------------------------------------------------------------------------
// helper macro for checking objects
//------------------------------------------------------------------------------

#define CHECK_MAGIC(object,kind)                                        \
{                                                                       \
    switch (object->magic)                                              \
    {                                                                   \
        case MAGIC :                                                    \
            /* the object is valid */                                   \
            break ;                                                     \
                                                                        \
        case FREED :                                                    \
            /* dangling pointer! */                                     \
            if (pr > 0) printf ("already freed!\n") ;                   \
            return (ERROR (GrB_UNINITIALIZED_OBJECT, (LOG,              \
                "%s is freed: [%s]", kind, name))) ;                    \
                                                                        \
        case MAGIC2 :                                                   \
            /* invalid */                                               \
            if (pr > 0) printf ("invalid\n") ;                          \
            return (ERROR (GrB_INVALID_OBJECT, (LOG,                    \
                "%s is invalid: [%s]", kind, name))) ;                  \
                                                                        \
        default :                                                       \
            /* uninitialized */                                         \
            if (pr > 0) printf ("uninititialized\n") ;                  \
            return (ERROR (GrB_UNINITIALIZED_OBJECT, (LOG,              \
                "%s is uninitialized: [%s]", kind, name))) ;            \
    }                                                                   \
}

//------------------------------------------------------------------------------
// type-specific macros
//------------------------------------------------------------------------------

// returns the largest possible value for a given type
#define PLUS_INF(type)            \
    _Generic                      \
    (                             \
        (type),                   \
        bool     : true         , \
        int8_t   : INT8_MAX     , \
        uint8_t  : UINT8_MAX    , \
        int16_t  : INT16_MAX    , \
        uint16_t : UINT16_MAX   , \
        int32_t  : INT32_MAX    , \
        uint32_t : UINT32_MAX   , \
        int64_t  : INT64_MAX    , \
        uint64_t : UINT64_MAX   , \
        float    : INFINITY     , \
        double   : INFINITY       \
    )

// returns the smallest possible value for a given type
#define MINUS_INF(type)           \
    _Generic                      \
    (                             \
        (type),                   \
        bool     : false        , \
        int8_t   : INT8_MIN     , \
        uint8_t  : 0            , \
        int16_t  : INT16_MIN    , \
        uint16_t : 0            , \
        int32_t  : INT32_MIN    , \
        uint32_t : 0            , \
        int64_t  : INT64_MIN    , \
        uint64_t : 0            , \
        float    : -INFINITY    , \
        double   : -INFINITY      \
    )

// true if the type is signed
#define IS_SIGNED(type)           \
    _Generic                      \
    (                             \
        (type),                   \
        bool     : false        , \
        int8_t   : true         , \
        uint8_t  : false        , \
        int16_t  : true         , \
        uint16_t : false        , \
        int32_t  : true         , \
        uint32_t : false        , \
        int64_t  : true         , \
        uint64_t : false        , \
        float    : true         , \
        double   : true           \
    )

// true if the type is integer (boolean is not integer)
#define IS_INTEGER(type)          \
    _Generic                      \
    (                             \
        (type),                   \
        bool     : false        , \
        int8_t   : true         , \
        uint8_t  : true         , \
        int16_t  : true         , \
        uint16_t : true         , \
        int32_t  : true         , \
        uint32_t : true         , \
        int64_t  : true         , \
        uint64_t : true         , \
        float    : false        , \
        double   : false          \
    )

// true if the type is boolean
#define IS_BOOLEAN(type)          \
    _Generic                      \
    (                             \
        (type),                   \
        bool     : true         , \
        int8_t   : false        , \
        uint8_t  : false        , \
        int16_t  : false        , \
        uint16_t : false        , \
        int32_t  : false        , \
        uint32_t : false        , \
        int64_t  : false        , \
        uint64_t : false        , \
        float    : false        , \
        double   : false          \
    )

// true if the type is float or double
#define IS_FLOAT(type)            \
    _Generic                      \
    (                             \
        (type),                   \
        bool     : false        , \
        int8_t   : false        , \
        uint8_t  : false        , \
        int16_t  : false        , \
        uint16_t : false        , \
        int32_t  : false        , \
        uint32_t : false        , \
        int64_t  : false        , \
        uint64_t : false        , \
        float    : true         , \
        double   : true           \
    )

//------------------------------------------------------------------------------
// division by zero
//------------------------------------------------------------------------------

// Integer division by zero is done the same way it's done in MATLAB.  This
// approach allows GraphBLAS to not terminate the user's application on
// divide-by-zero, and allows GraphBLAS results to be tested against MATLAB.
// To compute X/0: if X is zero, the result is zero (like NaN).  if X is
// negative, the result is the negative integer with biggest magnitude (like
// -infinity).  if X is positive, the result is the biggest positive integer
// (like +infinity).

// For places affected by this decision in the code do:
// grep "integer division"

// signed integer division x/zero: special cases for 0/0, -(neg)/0, (pos)/0
#define SIGNED_INTEGER_DIVISION_BY_ZERO(x)                                  \
(                                                                           \
    ((x) == 0) ?                                                            \
    (                                                                       \
        /* zero divided by zero gives 'Nan' */                              \
        0                                                                   \
    )                                                                       \
    :                                                                       \
    (                                                                       \
        /* x/0 and x is nonzero */                                          \
        ((x) < 0) ?                                                         \
        (                                                                   \
            /* x is negative: x/0 gives '-Inf' */                           \
            MINUS_INF (x)                                                   \
        )                                                                   \
        :                                                                   \
        (                                                                   \
            /* x is positive: x/0 gives '+Inf' */                           \
            PLUS_INF (x)                                                    \
        )                                                                   \
    )                                                                       \
)

// signed integer division x/y: special cases for y=-1 and y=0
#define SIGNED_INTEGER_DIVISION(x,y)                                        \
(                                                                           \
    ((y) == -1) ?                                                           \
    (                                                                       \
        /* INT32_MIN/(-1) causes floating point exception; avoid it  */     \
        -(x)                                                                \
    )                                                                       \
    :                                                                       \
    (                                                                       \
        ((y) == 0) ?                                                        \
        (                                                                   \
            /* x/0 */                                                       \
            SIGNED_INTEGER_DIVISION_BY_ZERO (x)                             \
        )                                                                   \
        :                                                                   \
        (                                                                   \
            /* normal case for signed integer division */                   \
            (x) / (y)                                                       \
        )                                                                   \
    )                                                                       \
)

// unsigned integer division x/zero: special cases for 0/0 and (pos)/0
#define UNSIGNED_INTEGER_DIVISION_BY_ZERO(x)                                \
(                                                                           \
    ((x) == 0) ?                                                            \
    (                                                                       \
        /* zero divided by zero gives 'Nan' */                              \
        0                                                                   \
    )                                                                       \
    :                                                                       \
    (                                                                       \
        /* x is positive: x/0 gives '+Inf' */                               \
        PLUS_INF (x)                                                        \
    )                                                                       \
)                                                                           \

// unsigned integer division x/y: special case for y=0
#define UNSIGNED_INTEGER_DIVISION(x,y)                                      \
(                                                                           \
    ((y) == 0) ?                                                            \
    (                                                                       \
        /* x/0 */                                                           \
        UNSIGNED_INTEGER_DIVISION_BY_ZERO (x)                               \
    )                                                                       \
    :                                                                       \
    (                                                                       \
        /* normal case for unsigned integer division */                     \
        (x) / (y)                                                           \
    )                                                                       \
)

// x/y when x and y are signed or unsigned integers
#define IDIV(x,y)                                                           \
(                                                                           \
    (IS_SIGNED (x)) ?                                                       \
    (                                                                       \
        SIGNED_INTEGER_DIVISION (x,y)                                       \
    )                                                                       \
    :                                                                       \
    (                                                                       \
        UNSIGNED_INTEGER_DIVISION (x,y)                                     \
    )                                                                       \
)

// 1/y when y is signed or unsigned integer
#define IMINV(y)                                                            \
(                                                                           \
    ((y) == 0) ?                                                            \
    (                                                                       \
        /* 1/0 */                                                           \
        PLUS_INF (y)                                                        \
    )                                                                       \
    :                                                                       \
    (                                                                       \
        /* normal case for integer minv, signed or unsigned */              \
        1 / (y)                                                             \
    )                                                                       \
)

// GraphBLAS includes a built-in GrB_DIV_BOOL operator, so boolean division
// must be defined.  There is no MATLAB equivalent since x/y for logical x and
// y is not permitted in MATLAB.  ANSI C11 does not provide a definition
// either, and dividing by zero (boolean false) will typically terminate an
// application.  In this GraphBLAS implementation, boolean division is treated
// as if it were int1, where 1/1 = 1, 0/1 = 0, 0/0 = integer NaN = 0, 1/0 =
// +infinity = 1.  Thus Z=X/Y is Z=X.  This is arbitrary, but it allows all
// operators to work on all types without causing run time exceptions.  It also
// means that GrB_DIV(x,y) is the same as GrB_FIRST(x,y) for boolean x and y.
// See for example GB_boolean_rename and Template/GB_ops_template.c.
// Similarly, GrB_MINV_BOOL, which is 1/x, is simply 'true' for all x.

//------------------------------------------------------------------------------
// typecasting
//------------------------------------------------------------------------------

// The ANSI C11 language specification states that results are undefined when
// typecasting a float or double to an integer value that is outside the range
// of the integer type.  However, most implementations provide a reasonable and
// repeatable result using modular arithmetic, just like conversions between
// integer types, so this is not changed here.  MATLAB uses a different
// strategy; any value outside the range is maxed out the largest or smallest
// integer.  Users of a C library would not expect this behavior.  They would
// expect instead whatever their C compiler would do in this case, even though
// it is technically undefined by the C11 standard.  Thus, this implementation
// of GraphBLAS does not attempt to second-guess the compiler when converting
// large floating-point values to integers.

// However, Inf's and NaN's are very unpredictable.  The same C compiler can
// generate different results with different optimization levels.  Only bool is
// defined by the ANSI C11 standard (NaN converts to true, since Nan != 0 is
// true).

// This unpredictability with Inf's and NaN's causes the GraphBLAS tests to
// fail in unpredictable ways.  Therefore, in this implementation of GraphBLAS,
// a float or double +Inf is converted to the largest integer, -Inf to the
// smallest integer, and NaN to zero.  This is the same behavior as MATLAB.

// cast a value x to the type of z
#define CAST(z,x)                                                           \
{                                                                           \
    if (IS_INTEGER (z) && IS_FLOAT (x))                                     \
    {                                                                       \
        switch (fpclassify ((double) (x)))                                  \
        {                                                                   \
            case FP_ZERO:                                                   \
            case FP_NORMAL:                                                 \
            case FP_SUBNORMAL:                                              \
                z = (x) ;                                                   \
                break ;                                                     \
            case FP_NAN:                                                    \
                z = 0 ;                                                     \
                break ;                                                     \
            case FP_INFINITE:                                               \
                z = ((x) > 0) ? PLUS_INF (z) : MINUS_INF (z) ;              \
                break ;                                                     \
        }                                                                   \
    }                                                                       \
    else                                                                    \
    {                                                                       \
        /* trust the built-in typecast */                                   \
        z = (x) ;                                                           \
    }                                                                       \
}


//------------------------------------------------------------------------------
// GB_BINARY_SEARCH
//------------------------------------------------------------------------------

// search for integer i in the list X [pleft...pright]; no zombies.
// The list X [pleft ... pright] is in ascending order.  It may have
// duplicates.

#define GB_BINARY_TRIM_SEARCH(i,X,pleft,pright)                             \
{                                                                           \
    /* binary search of X [pleft ... pright] for integer i */               \
    while (pleft < pright)                                                  \
    {                                                                       \
        int64_t pmiddle = (pleft + pright) / 2 ;                            \
        if (X [pmiddle] < i)                                                \
        {                                                                   \
            /* if in the list, it appears in [pmiddle+1..pright] */         \
            pleft = pmiddle + 1 ;                                           \
        }                                                                   \
        else                                                                \
        {                                                                   \
            /* if in the list, it appears in [pleft..pmiddle] */            \
            pright = pmiddle ;                                              \
        }                                                                   \
    }                                                                       \
    /* binary search is narrowed down to a single item */                   \
    /* or it has found the list is empty */                                 \
    ASSERT (pleft == pright || pleft == pright + 1) ;                       \
}

// GB_BINARY_SEARCH:
// If found is true then X [pleft == pright] == i.  If duplicates appear then
// X [pleft] is any one of the entries with value i in the list.
// If found is false then
//    X [original_pleft ... pleft-1] < i and
//    X [pleft+1 ... original_pright] > i holds.
// The value X [pleft] may be either < or > i.
#define GB_BINARY_SEARCH(i,X,pleft,pright,found)                            \
{                                                                           \
    GB_BINARY_TRIM_SEARCH (i, X, pleft, pright) ;                           \
    found = (pleft == pright && X [pleft] == i) ;                           \
}

// GB_BINARY_SPLIT_SEARCH
// If found is true then X [pleft] == i.  If duplicates appear then X [pleft]
//    is any one of the entries with value i in the list.
// If found is false then
//    X [original_pleft ... pleft-1] < i and
//    X [pleft ... original_pright-1] > i holds, and pleft-1 == pright
// If X has no duplicates, then whether or not i is found,
//    X [original_pleft ... pleft-1] < i and
//    X [pleft ... original_pright-1] >= i holds.
#define GB_BINARY_SPLIT_SEARCH(i,X,pleft,pright,found)                      \
{                                                                           \
    GB_BINARY_SEARCH (i, X, pleft, pright, found)                           \
    if (!found && (pleft == pright))                                        \
    {                                                                       \
        if (i > X [pleft])                                                  \
        {                                                                   \
            pleft++ ;                                                       \
        }                                                                   \
        else                                                                \
        {                                                                   \
            pright++ ;                                                      \
        }                                                                   \
    }                                                                       \
}

//------------------------------------------------------------------------------
// GB_BINARY_ZOMBIE
//------------------------------------------------------------------------------

#define GB_BINARY_TRIM_ZOMBIE(i,X,pleft,pright)                             \
{                                                                           \
    /* binary search of X [pleft ... pright] for integer i */               \
    while (pleft < pright)                                                  \
    {                                                                       \
        int64_t pmiddle = (pleft + pright) / 2 ;                            \
        if (i > UNFLIP (X [pmiddle]))                                       \
        {                                                                   \
            /* if in the list, it appears in [pmiddle+1..pright] */         \
            pleft = pmiddle + 1 ;                                           \
        }                                                                   \
        else                                                                \
        {                                                                   \
            /* if in the list, it appears in [pleft..pmiddle] */            \
            pright = pmiddle ;                                              \
        }                                                                   \
    }                                                                       \
    /* binary search is narrowed down to a single item */                   \
    /* or it has found the list is empty */                                 \
    ASSERT (pleft == pright || pleft == pright + 1) ;                       \
}

#define GB_BINARY_ZOMBIE(i,X,pleft,pright,found,nzombies,is_zombie)         \
{                                                                           \
    if (nzombies > 0)                                                       \
    {                                                                       \
        GB_BINARY_TRIM_ZOMBIE (i, X, pleft, pright) ;                       \
        found = false ;                                                     \
        is_zombie = false ;                                                 \
        if (pleft == pright)                                                \
        {                                                                   \
            int64_t i2 = X [pleft] ;                                        \
            is_zombie = IS_ZOMBIE (i2) ;                                    \
            if (is_zombie)                                                  \
            {                                                               \
                i2 = FLIP (i2) ;                                            \
            }                                                               \
            found = (i == i2) ;                                             \
        }                                                                   \
    }                                                                       \
    else                                                                    \
    {                                                                       \
        is_zombie = false ;                                                 \
        GB_BINARY_SEARCH(i,X,pleft,pright,found)                            \
    }                                                                       \
}

//------------------------------------------------------------------------------
// index lists I and J
//------------------------------------------------------------------------------

// kind of index list, Ikind and Jkind:
#define GB_ALL 0
#define GB_RANGE 1
#define GB_STRIDE 2
#define GB_LIST 4

void GB_ijlength            // get the length and kind of an index list I
(
    const GrB_Index *I,     // list of indices (actual or implicit)
    const int64_t ni,       // length I, or special
    const int64_t limit,    // indices must be in the range 0 to limit-1
    int64_t *nI,            // actual length of I
    int *Ikind,             // kind of I: GB_ALL, GB_RANGE, GB_STRIDE, GB_LIST
    int64_t Icolon [3]      // begin:inc:end for all but GB_LIST
) ;

GrB_Info GB_ijproperties        // check I and determine its properties
(
    const GrB_Index *I,         // list of indices, or GrB_ALL
    const int64_t ni,           // length I, or special
    const int64_t nI,           // actual length from GB_ijlength
    const int64_t limit,        // I must be in the range 0 to limit-1
    const int64_t Ikind,        // kind of I, from GB_ijlength
    const int64_t Icolon [3],   // begin:inc:end from GB_ijlength
    bool *I_is_unsorted,        // true if I is out of order
    bool *I_is_contig,          // true if I is a contiguous list, imin:imax
    int64_t *imin_result,       // min (I)
    int64_t *imax_result,       // max (I)
    bool is_I
) ;

GrB_Info GB_ijsort
(
    const GrB_Index *I, // index array of size ni
    int64_t *p_ni,      // input: size of I, output: number of indices in I2
    GrB_Index **p_I2    // output array of size ni, where I2 [0..ni2-1]
                        // contains the sorted indices with duplicates removed.
) ;

// given k, return the kth item i = I [k] in the list
static inline int64_t GB_ijlist     // get the kth item in a list of indices
(
    const GrB_Index *I,         // list of indices
    const int64_t k,            // return i = I [k], the kth item in the list
    const int Ikind,            // GB_ALL, GB_RANGE, GB_STRIDE, or GB_LIST
    const int64_t Icolon [3]    // begin:inc:end for all but GB_LIST
)
{
    if (Ikind == GB_ALL)
    { 
        // I is ":"
        return (k) ;
    }
    else if (Ikind == GB_RANGE)
    { 
        // I is begin:end
        return (Icolon [GxB_BEGIN] + k) ;
    }
    else if (Ikind == GB_STRIDE)
    { 
        // I is begin:inc:end
        // note that iinc can be negative or even zero
        return (Icolon [GxB_BEGIN] + k * Icolon [GxB_INC]) ;
    }
    else // Ikind == GB_LIST
    { 
        ASSERT (Ikind == GB_LIST) ;
        ASSERT (I != NULL) ;
        return (I [k]) ;
    }
}

// given i and I, return true there is a k so that i is the kth item in I
static inline bool GB_ij_is_in_list // determine if i is in the list I
(
    const GrB_Index *I,         // list of indices for GB_LIST
    const int64_t nI,           // length of I
    int64_t i,                  // find i = I [k] in the list
    const int Ikind,            // GB_ALL, GB_RANGE, GB_STRIDE, or GB_LIST
    const int64_t Icolon [3]    // begin:inc:end for all but GB_LIST
)
{
    if (Ikind == GB_ALL)
    { 
        // I is ":", all indices are in the sequence
        return (true) ;
    }
    else if (Ikind == GB_RANGE)
    { 
        // I is begin:end
        int64_t b = Icolon [GxB_BEGIN] ;
        int64_t e = Icolon [GxB_END] ;
        if (i < b) return (false) ;
        if (i > e) return (false) ;
        return (true) ;
    }
    else if (Ikind == GB_STRIDE)
    {
        // I is begin:inc:end
        // note that inc can be negative or even zero
        int64_t b   = Icolon [GxB_BEGIN] ;
        int64_t inc = Icolon [GxB_INC] ;
        int64_t e   = Icolon [GxB_END] ;
        if (inc == 0)
        { 
            // I is empty if inc is zero, so i is not in I
            return (false) ;
        }
        else if (inc > 0)
        { 
            // forward direction, increment is positive
            // I = b:inc:e = [b, b+inc, b+2*inc, ..., e]
            if (i < b) return (false) ;
            if (i > e) return (false) ;
            // now i is in the range [b ... e]
            ASSERT (b <= i && i <= e) ;
            i = i - b ;
            ASSERT (0 <= i && i <= (e-b)) ;
            // the sequence I-b = [0, inc, 2*inc, ... e-b].
            // i is in the sequence if i % inc == 0
            return (i % inc == 0) ;
        }
        else // inc < 0
        { 
            // backwards direction, increment is negative
            inc = -inc ;
            // now inc is positive
            ASSERT (inc > 0) ;
            // I = b:(-inc):e = [b, b-inc, b-2*inc, ... e]
            if (i > b) return (false) ;
            if (i < e) return (false) ;
            // now i is in the range of the sequence, [b down to e]
            ASSERT (e <= i && i <= b) ;
            i = b - i ;
            ASSERT (0 <= i && i <= (b-e)) ;
            // b-I = 0:(inc):(b-e) = [0, inc, 2*inc, ... (b-e)]
            // i is in the sequence if i % inc == 0
            return (i % inc == 0) ;
        }
    }
    else // Ikind == GB_LIST
    { 
        ASSERT (Ikind == GB_LIST) ;
        ASSERT (I != NULL) ;
        // search for i in the sorted list I
        bool found ;
        int64_t pleft = 0 ;
        int64_t pright = nI-1 ;
        GB_BINARY_SEARCH (i, I, pleft, pright, found) ;
        return (found) ;
    }
}

//------------------------------------------------------------------------------
// GB_bracket_left
//------------------------------------------------------------------------------

// Given a sorted list X [kleft:kright], and a range imin:..., modify kleft so
// that the smaller sublist X [kleft:kright] contains the range imin:...

static inline void GB_bracket_left
(
    const int64_t imin,
    const int64_t *restrict X,  // input list is in X [kleft:kright]
    int64_t *kleft,
    const int64_t kright
)
{
    // tighten kleft
    int64_t len = kright - (*kleft) + 1 ;
    if (len > 0 && X [(*kleft)] < imin)
    { 
        // search for imin in X [kleft:kright]
        int64_t pleft = (*kleft) ;
        int64_t pright = kright ;
        GB_BINARY_TRIM_SEARCH (imin, X, pleft, pright) ;
        (*kleft) = pleft ;
    }
}

//------------------------------------------------------------------------------
// GB_bracket_right
//------------------------------------------------------------------------------

// Given a sorted list X [kleft:kright], and a range ...:imax, modify kright so
// that the smaller sublist X [kleft:kright] contains the range ...:imax.

static inline void GB_bracket_right
(
    const int64_t imax,
    const int64_t *restrict X,  // input list is in X [kleft:kright]
    const int64_t kleft,
    int64_t *kright
)
{
    // tighten kright
    int64_t len = (*kright) - kleft + 1 ;
    if (len > 0 && imax < X [(*kright)])
    { 
        // search for imax in X [kleft:kright]
        int64_t pleft = kleft ;
        int64_t pright = (*kright) ;
        GB_BINARY_TRIM_SEARCH (imax, X, pleft, pright) ;
        (*kright) = pleft ;
    }
}

//------------------------------------------------------------------------------
// GB_bracket
//------------------------------------------------------------------------------

// Given a sorted list X [kleft:kright], and a range imin:imax, find the
// kleft_new and kright_new so that the smaller sublist X
// [kleft_new:kright_new] contains the range imin:imax.

// Zombies are not tolerated.

static inline void GB_bracket
(
    const int64_t imin,         // search for entries in the range imin:imax
    const int64_t imax,
    const int64_t *restrict X,  // input list is in X [kleft:kright]
    const int64_t kleft_in,
    const int64_t kright_in,
    int64_t *kleft_new,         // output list is in X [kleft_new:kright_new]
    int64_t *kright_new
)
{ 

    int64_t kleft  = kleft_in ;
    int64_t kright = kright_in ;

    if (imin > imax)
    { 
        // imin:imax is empty, make X [kleft:kright] empty
        (*kleft_new ) = kleft ;
        (*kright_new) = kleft-1 ;
        return ;
    }

    // Find kleft and kright so that X [kleft:kright] contains all of imin:imax

    // tighten kleft
    GB_bracket_left (imin, X, &kleft, kright) ;

    // tighten kright
    GB_bracket_right (imax, X, kleft, &kright) ;

    // list has been trimmed
    ASSERT (IMPLIES (kleft > kleft_in, X [kleft-1] < imin)) ;
    ASSERT (IMPLIES (kright < kright_in, imax < X [kright+1])) ;

    // return result
    (*kleft_new ) = kleft ;
    (*kright_new) = kright ;
}

//------------------------------------------------------------------------------
// NaN handling
//------------------------------------------------------------------------------

// For floating-point computations, SuiteSparse:GraphBLAS relies on the IEEE
// 754 standard for the basic operations (+ - / *).  Comparison operators also
// work as they should; any comparison with NaN is always false, even
// eq(NaN,NaN) is false.  This follows the IEEE 754 standard.

// For integer MIN and MAX, SuiteSparse:GraphBLAS relies on one comparison:

// z = min(x,y) = (x < y) ? x : y
// z = max(x,y) = (x > y) ? x : y

// However, this is not suitable for floating-point x and y.  Comparisons with
// NaN always return false, so if either x or y are NaN, then z = y, for both
// min(x,y) and max(x,y).  In MATLAB, min(3,NaN), min(NaN,3), max(3,NaN), and
// max(NaN,3) are all 3, which is another interpretation.  The MATLAB min and
// max functions have a 3rd argument that specifies how NaNs are handled:
// 'omitnan' (default) and 'includenan'.  In SuiteSparse:GraphBLAS, the min and
// max functions are the same as 'includenan' in MATLAB.

// Below is a complete comparison of MATLAB and GraphBLAS.  Both tables are the
// results for both min and max (they return the same results in these cases):

//   x    y  MATLAB    MATLAB   (x<y)?x:y   SuiteSparse:
//           omitnan includenan             GraphBLAS
//   3    3     3        3          3        3
//   3   NaN    3       NaN        NaN      NaN
//  NaN   3     3       NaN         3       NaN
//  NaN  NaN   NaN      NaN        NaN      NaN

// The 'includenan' result is the one SuiteSparse:GraphBLAS relies on.  The
// result min(x,y) and max(x,y) should be NaN if either x or y are NaN.  NaN
// means 'not any number', and any computations should propagate the result of
// NaN if a function depends on that value.  If x = 5 then z=min(3,x) is 3, but
// if x = 2, z=min(3,x) is 2.  Since z depends on x, z should be NaN if x is
// NaN.

// suitable for all types, but given unique names to match IMIN, FMIN:
#define IABS(x) (((x) >= 0) ? (x) : (-(x)))
#define FABS(x) (((x) >= 0) ? (x) : (-(x)))

// suitable for integers, and non-NaN floating point:
#define IMAX(x,y) (((x) > (y)) ? (x) : (y))
#define IMIN(x,y) (((x) < (y)) ? (x) : (y))

// MIN for floating-point, same as min(x,y,'includenan') in MATLAB
#define FMIN(x,y) ((isnan (x) || isnan (y)) ? NAN : IMIN (x,y))

// MAX for floating-point, same as max(x,y,'includenan') in MATLAB
#define FMAX(x,y) ((isnan (x) || isnan (y)) ? NAN : IMAX (x,y))

//==============================================================================
// GrB_Matrix iterator and related functions
//==============================================================================

// find k so that j == Ah [k], if it appears in the list.  k is not needed
// by the caller, just the variables pstart, pend, pleft, and found.

static inline bool GB_lookup        // find j = Ah [k] in the hyperlist of a
(
    const bool A_is_hyper,          // true if A is hypersparse
    const int64_t *restrict Ah,     // A->h [0..A->nvec-1]: list of vectors
    const int64_t *restrict Ap,     // A->p [0..A->nvec  ]: pointers to vectors
    int64_t *restrict pleft,        // look only in A->h [pleft..pright]
    int64_t pright,                 // normally A->nvec-1, but can be trimmed
//  const int64_t nvec,             // A->nvec: number of vectors
    const int64_t j,                // vector to find, as j = Ah [k]
    int64_t *restrict pstart,       // start of vector: Ap [k]
    int64_t *restrict pend          // end of vector: Ap [k+1]
)
{
    if (A_is_hyper)
    {
        // to search the whole Ah list, use on input:
        // pleft = 0 ; pright = nvec-1 ;
        bool found ;
        GB_BINARY_SEARCH (j, Ah, (*pleft), pright, found) ;
        if (found)
        { 
            // j appears in the hyperlist at Ah [pleft]
            // k = (*pleft)
            (*pstart) = Ap [(*pleft)] ;
            (*pend)   = Ap [(*pleft)+1] ;
        }
        else
        { 
            // j does not appear in the hyperlist Ah
            // k = -1
            (*pstart) = -1 ;
            (*pend)   = -1 ;
        }
        return (found) ;
    }
    else
    { 
        // A is not hypersparse; j always appears
        // k = j
        (*pstart) = Ap [j] ;
        (*pend)   = Ap [j+1] ;
        return (true) ;
    }
}

//------------------------------------------------------------------------------
// GBI_iterator: iterate over the vectors one, two, or three matrices
//------------------------------------------------------------------------------

// Any matrix of the three matrices may be standard sparse or hypersparse.  The
// A->h and A->p components of the matrices may not change during the
// iteration.

// The Iter->* content of a GBI_iterator Iter is accessed only in this file.

// All typedefs, functions, and macros that operate on the
// SuiteSparse:GraphBLAS iterator have names that start with the GBI prefix.

typedef struct
{
    bool any_hyper ;            // true if any matrix is hypersparse
    bool is_hyper [3] ;         // true if A [0..2] is hypersparse

    int64_t j ;                 // name of current vector
    int64_t vdim ;              // A->vdim: number of vectors in A [0..2]

    const int64_t *p [3] ;      // vector pointers A->p of A [0..2]
    int64_t pstart [3] ;        // start of current vector
    int64_t pend [3] ;          // one past the end of current vector

    int64_t nvec [3] ;          // A->nvec: number of vectors in A [0..2]
    const int64_t *h [3] ;      // A->h: hyperlist of vectors in A [0..2]
    int64_t k [3] ;             // current index into hyperlist h [0..2]
    int64_t jj [3] ;            // current vector jj [.] = h [k [.]]

} GBI_iterator ;

//------------------------------------------------------------------------------
// GBI#_init
//------------------------------------------------------------------------------

// Initialize an iterator for one, two, or three matrices.  Extracts the vector
// content (is_hyper, vdim, p, h, and nvec) of A, B, and/or C and caches them
// in the GBI_iterator Iter, and initializes the counters Iter->j and Iter->k
// [0..2].  The GBI#_init functions are the only inline functions that directly
// access the matrix content A->(...) for A, B, and/or C.  All other iterator
// functions and macros access the cached copies in Iter.

static inline void GBI3_init
(
    GBI_iterator *Iter,
    const GrB_Matrix A,
    const GrB_Matrix B,
    const GrB_Matrix C
)
{ 
    // A, B, and C must have the same vdim.  vlen can vary
    ASSERT (A->vdim == B->vdim) ;
    ASSERT (A->vdim == C->vdim) ;

    Iter->any_hyper = (A->is_hyper || B->is_hyper || C->is_hyper) ;
    Iter->j = 0 ;
    Iter->vdim = A->vdim ;

    Iter->p [0] = A->p ;
    Iter->p [1] = B->p ;
    Iter->p [2] = C->p ;

    Iter->h [0] = A->h ;
    Iter->h [1] = B->h ;
    Iter->h [2] = C->h ;

    Iter->jj [0] = 0 ;
    Iter->jj [1] = 0 ;
    Iter->jj [2] = 0 ;

    Iter->k [0] = 0 ;
    Iter->k [1] = 0 ;
    Iter->k [2] = 0 ;

    Iter->pstart [0] = 0 ;
    Iter->pstart [1] = 0 ;
    Iter->pstart [2] = 0 ;

    Iter->pend [0] = 0 ;
    Iter->pend [1] = 0 ;
    Iter->pend [2] = 0 ;

    Iter->nvec [0] = A->nvec ;
    Iter->nvec [1] = B->nvec ;
    Iter->nvec [2] = C->nvec ;

    Iter->is_hyper [0] = A->is_hyper ;
    Iter->is_hyper [1] = B->is_hyper ;
    Iter->is_hyper [2] = C->is_hyper ;

}

static inline void GBI2_init
(
    GBI_iterator *Iter,
    const GrB_Matrix A,
    const GrB_Matrix B
)
{ 
    // A and B must have the same vdim.  vlen can vary
    ASSERT (A->vdim == B->vdim) ;

    Iter->any_hyper = (A->is_hyper || B->is_hyper) ;
    Iter->j = 0 ;
    Iter->vdim = A->vdim ;
    Iter->p [0] = A->p ;
    Iter->p [1] = B->p ;
//  Iter->p [2] = NULL ;

    Iter->h [0] = A->h ;
    Iter->h [1] = B->h ;
//  Iter->h [2] = NULL ;

    Iter->jj [0] = 0 ;
    Iter->jj [1] = 0 ;
//  Iter->jj [2] = 0 ;

    Iter->k [0] = 0 ;
    Iter->k [1] = 0 ;
//  Iter->k [2] = 0 ;

    Iter->pstart [0] = 0 ;
    Iter->pstart [1] = 0 ;
//  Iter->pstart [2] = 0 ;

    Iter->pend [0] = 0 ;
    Iter->pend [1] = 0 ;
//  Iter->pend [2] = 0 ;

    Iter->nvec [0] = A->nvec ;
    Iter->nvec [1] = B->nvec ;
//  Iter->nvec [2] = 0 ;

    Iter->is_hyper [0] = A->is_hyper ;
    Iter->is_hyper [1] = B->is_hyper ;
//  Iter->is_hyper [2] = false ;

}

static inline void GBI1_init
(
    GBI_iterator *Iter,
    const GrB_Matrix A
)
{ 

    Iter->any_hyper = (A->is_hyper) ;
    Iter->j = 0 ;
    Iter->vdim = A->vdim ;
    Iter->p [0] = A->p ;
//  Iter->p [1] = NULL ;
//  Iter->p [2] = NULL ;

    Iter->h [0] = A->h ;
//  Iter->h [1] = NULL ;
//  Iter->h [2] = NULL ;

    Iter->jj [0] = 0 ;
//  Iter->jj [1] = 0 ;
//  Iter->jj [2] = 0 ;

    Iter->k [0] = 0 ;
//  Iter->k [1] = 0 ;
//  Iter->k [2] = 0 ;

    Iter->pstart [0] = 0 ;
//  Iter->pstart [1] = 0 ;
//  Iter->pstart [2] = 0 ;

    Iter->pend [0] = 0 ;
//  Iter->pend [1] = 0 ;
//  Iter->pend [2] = 0 ;

    Iter->nvec [0] = A->nvec ;
//  Iter->nvec [1] = 0 ;
//  Iter->nvec [2] = 0 ;

    Iter->is_hyper [0] = A->is_hyper ;
//  Iter->is_hyper [1] = false ;
//  Iter->is_hyper [2] = false ;

}

static inline void GBI3s_init
(
    GBI_iterator *Iter,
    const GrB_Matrix A,
    const GrB_Matrix B
)
{ 
    GBI2_init (Iter, A, B) ;
}

static inline void GBI2s_init
(
    GBI_iterator *Iter,
    const GrB_Matrix A
)
{ 
    GBI1_init (Iter, A) ;
}

//------------------------------------------------------------------------------
// GBI#_while
//------------------------------------------------------------------------------

// Test the condition for a loop over one, two, or three matrices and
// prepare to access the vector at each matrix.  The inline
// GBI#_while functions are only used in the each_vector* macros.

// get the index of the next vector in a given matrix
#define GBI_jjget(Iter,matrix)                                              \
{                                                                           \
    if (!(Iter->is_hyper [matrix]))                                         \
    {                                                                       \
        /* the matrix is not hypersparse */                                 \
        Iter->jj [matrix] = Iter->k [matrix] ;                              \
    }                                                                       \
    else if (Iter->k [matrix] < Iter->nvec [matrix])                        \
    {                                                                       \
        /* jj is the kth non-empty vector in the matrix */                  \
        Iter->jj [matrix] = (Iter->h [matrix]) [Iter->k [matrix]] ;         \
    }                                                                       \
    else                                                                    \
    {                                                                       \
        /* list of non-empty vectors exhausted; return sentinel value */    \
        Iter->jj [matrix] = Iter->vdim ;                                    \
    }                                                                       \
}

// get the pointers to the next vector A (:,j), not hypersparse
#define GBI_pget(Iter,matrix)                                               \
{                                                                           \
    Iter->pstart [matrix] = (Iter->p [matrix]) [Iter->j  ] ;                \
    Iter->pend   [matrix] = (Iter->p [matrix]) [Iter->j+1] ;                \
}

// get the pointers to the kth vector A (:,j), hypersparse
#define GBI_hget(Iter,matrix)                                               \
{                                                                           \
    Iter->pstart [matrix] = (Iter->p [matrix]) [Iter->k [matrix]  ] ;       \
    Iter->pend   [matrix] = (Iter->p [matrix]) [Iter->k [matrix]+1] ;       \
}

// Get the pointers to the next vector A (:,j), sparse or hypersparse.  The
// vector appears in Ai and Ax [Iter->pstart [matrix]...Iter->pend [matrix]-1].
#define GBI_phget(Iter,matrix)                                              \
{                                                                           \
    if (!(Iter->is_hyper [matrix]))                                         \
    {                                                                       \
        /* A (:,j) is the jth vector in the non-hypersparse matrix */       \
        GBI_pget (Iter, matrix) ;                                           \
    }                                                                       \
    else if (Iter->j == Iter->jj [matrix])                                  \
    {                                                                       \
        /* A (:,j) is the kth vector in the hypersparse matrix */           \
        GBI_hget (Iter, matrix) ;                                           \
    }                                                                       \
    else                                                                    \
    {                                                                       \
        /* A (:,j) does not appear in the hypersparse matrix */             \
        Iter->pstart [matrix] = -1 ;                                        \
        Iter->pend   [matrix] = -1 ;                                        \
    }                                                                       \
}

static inline bool GBI3_while
(
    GBI_iterator *Iter
)
{

    if (Iter->any_hyper)
    { 
        // get next vector from A, B and/or C
        GBI_jjget (Iter, 0) ;
        GBI_jjget (Iter, 1) ;
        GBI_jjget (Iter, 2) ;
        Iter->j = IMIN (Iter->jj [0], Iter->jj [1]) ;
        Iter->j = IMIN (Iter->j     , Iter->jj [2]) ;

        // test if the end of the matrix has been reached
        if (Iter->j >= Iter->vdim) return (false) ;

        // get the start and end of each vector of A, B, and C
        GBI_phget (Iter, 0) ;
        GBI_phget (Iter, 1) ;
        GBI_phget (Iter, 2) ;
    }
    else
    { 
        // test if the end of the matrix has been reached
        if (Iter->j >= Iter->vdim) return (false) ;

        // get the start and end of the jth vector of A, B, and C
        GBI_pget (Iter, 0) ;
        GBI_pget (Iter, 1) ;
        GBI_pget (Iter, 2) ;
    }

    // end of matrix has not yet been reached
    return (true) ;
}

static inline bool GBI3s_while
(
    GBI_iterator *Iter
)
{

    // test if the end of the matrix has been reached
    if (Iter->j >= Iter->vdim) return (false) ;

    if (Iter->any_hyper)
    { 
        // get next vector from A and/or B
        GBI_jjget (Iter, 0) ;
        GBI_jjget (Iter, 1) ;
        GBI_phget (Iter, 0) ;
        GBI_phget (Iter, 1) ;
    }
    else
    { 
        // get the start and end of the jth vector of A and B
        GBI_pget (Iter, 0) ;
        GBI_pget (Iter, 1) ;
    }

    // end of matrix has not yet been reached
    return (true) ;
}

static inline bool GBI2_while
(
    GBI_iterator *Iter
)
{

    if (Iter->any_hyper)
    { 
        // get next vector from A and/or B
        GBI_jjget (Iter, 0) ;
        GBI_jjget (Iter, 1) ;
        Iter->j = IMIN (Iter->jj [0], Iter->jj [1]) ;

        // test if the end of the matrix has been reached
        if (Iter->j >= Iter->vdim) return (false) ;

        // get the start and end of each vector of A and B
        GBI_phget (Iter, 0) ;
        GBI_phget (Iter, 1) ;
    }
    else
    { 
        // test if the end of the matrix has been reached
        if (Iter->j >= Iter->vdim) return (false) ;

        // get the start and end of the jth vector of A and B
        GBI_pget (Iter, 0) ;
        GBI_pget (Iter, 1) ;
    }

    // end of matrix has not yet been reached
    return (true) ;
}

static inline bool GBI2s_while
(
    GBI_iterator *Iter
)
{

    // test if the end of the matrix has been reached
    if (Iter->j >= Iter->vdim) return (false) ;

    if (Iter->any_hyper)
    { 
        // get next vector from A
        GBI_jjget (Iter, 0) ;
        GBI_phget (Iter, 0) ;
    }
    else
    { 
        // get the start and end of the jth vector of A and B
        GBI_pget (Iter, 0) ;
    }

    // end of matrix has not yet been reached
    return (true) ;
}

static inline bool GBI1_while
(
    GBI_iterator *Iter
)
{

    if (Iter->any_hyper)
    { 
        // test if the end of the matrix has been reached
        if (Iter->k [0] >= Iter->nvec [0]) return (false) ;

        // get next vector from A
        Iter->j = (Iter->h [0]) [Iter->k [0]] ;

        // get the start and end of the next vector of A
        GBI_hget (Iter, 0) ;

    }
    else
    { 
        // test if the end of the matrix has been reached
        if (Iter->j >= Iter->vdim) return (false) ;

        // get the start and end of the jth vector of A
        GBI_pget (Iter, 0) ;
    }

    // end of matrix has not yet been reached
    return (true) ;
}

// The GBI_jjget, GBI_pget, GBI_hget, and GBI_phget macros
// are only directly used in the GBI#_while functions defined above.
#undef GBI_jjget
#undef GBI_pget
#undef GBI_hget
#undef GBI_phget

//------------------------------------------------------------------------------
// GBI#_next
//------------------------------------------------------------------------------

// Advance to the next vector of one, two, or three matrices.  The inline
// GBI#_next functions are only used in the each_vector* macros.

#define GBI_jnext(Iter,matrix)                                                \
{                                                                             \
    if (!(Iter->is_hyper [matrix]))                                           \
    {                                                                         \
        /* matrix is not hypersparse, so always advance to next vector */     \
        Iter->k [matrix]++ ;                                                  \
    }                                                                         \
    else if (Iter->j == Iter->jj [matrix])                                    \
    {                                                                         \
        /* matrix is hypersparse, and Ah [k] == j; advance to next vector */  \
        ASSERT (Iter->j == (Iter->h [matrix])[Iter->k [matrix]]) ;            \
        Iter->k [matrix]++ ;                                                  \
    }                                                                         \
}

static inline void GBI3_next
(
    GBI_iterator *Iter
)
{
    if (Iter->any_hyper)
    { 
        // advance to the next vector of A, B, and/or C
        GBI_jnext (Iter, 0) ;
        GBI_jnext (Iter, 1) ;
        GBI_jnext (Iter, 2) ;
    }
    else
    { 
        // A, B, and C are non-hypersparse
        Iter->j++ ;
    }
}

static inline void GBI3s_next
(
    GBI_iterator *Iter
)
{
    if (Iter->any_hyper)
    { 
        // advance to the next vector of A, and/or B
        GBI_jnext (Iter, 0) ;
        GBI_jnext (Iter, 1) ;
    }
    // C is an implicit scalar
    Iter->j++ ;
}

static inline void GBI2_next
(
    GBI_iterator *Iter
)
{
    if (Iter->any_hyper)
    { 
        // advance to the next vector of A and/or B
        GBI_jnext (Iter, 0) ;
        GBI_jnext (Iter, 1) ;
    }
    else
    { 
        // both A and B are non-hypersparse
        Iter->j++ ;
    }
}

static inline void GBI2s_next
(
    GBI_iterator *Iter
)
{
    if (Iter->any_hyper)
    { 
        // advance to the next vector of A
        GBI_jnext (Iter, 0) ;
    }
    // B is an implicit scalar
    Iter->j++ ;
}

static inline void GBI1_next
(
    GBI_iterator *Iter
)
{
    if (Iter->any_hyper)
    { 
        // always advance to the next vector of A
        Iter->k [0]++ ;
    }
    else
    { 
        // A is non-hypersparse
        Iter->j++ ;
    }
}

// the GBI_jnext* macros are only used in the GBI#_next functions above.
#undef GBI_jnext

//------------------------------------------------------------------------------
// GBI#_initj
//------------------------------------------------------------------------------

// Get the column at the current iteration, and the start/end pointers
// of column j in the three matrices A, B, and C.  These macros are used
// in the for (each_entry (...)) loop, and also at the start of the
// body of the for (each_vector (...)) loops, outside of this file.

#define GBI3_initj(Iter,j_,pstart0,pend0,pstart1,pend1,pstart2,pend2)  \
    j_      = Iter.j,                                                  \
    pstart0 = Iter.pstart [0],                                         \
    pend0   = Iter.pend   [0],                                         \
    pstart1 = Iter.pstart [1],                                         \
    pend1   = Iter.pend   [1],                                         \
    pstart2 = Iter.pstart [2],                                         \
    pend2   = Iter.pend   [2]

// get the column at the current iteration, and the start/end pointers
// of column j in the two matrices A and B

#define GBI2_initj(Iter,j_,pstart0,pend0,pstart1,pend1)                \
    j_      = Iter.j,                                                  \
    pstart0 = Iter.pstart [0],                                         \
    pend0   = Iter.pend   [0],                                         \
    pstart1 = Iter.pstart [1],                                         \
    pend1   = Iter.pend   [1]

// get the column at the current iteration, and the start/end pointers
// of column j in the matrix A

#define GBI1_initj(Iter,j_,pstart0,pend0)                              \
    j_      = Iter.j,                                                  \
    pstart0 = Iter.pstart [0],                                         \
    pend0   = Iter.pend   [0]

//------------------------------------------------------------------------------
// for-loop control: iterate over the vectors and entries of 1, 2, or 3 matrices
//------------------------------------------------------------------------------

// iterate over one matrix A
#define each_vector(Iter,A)                                                 \
    GBI1_init (&Iter,A)     ; GBI1_while (&Iter)  ; GBI1_next (&Iter)

#define for_each_vector(A)                                                  \
    GBI_iterator Iter ;                                                     \
    for (each_vector (Iter, A))

// iterate the union of two matrices A and B
#define each_vector2(Iter,A,B)                                              \
    GBI2_init (&Iter,A,B)   ; GBI2_while (&Iter)  ; GBI2_next (&Iter)

#define for_each_vector2(A,B)                                               \
    GBI_iterator Iter ;                                                     \
    for (each_vector2 (Iter, A, B))

// iterate the union of one or two matrices A and B.  B may be NULL.
#define each_vector2m(Iter,A,B)                                             \
    (B == NULL) ? GBI1_init  (&Iter,A) : GBI2_init  (&Iter,A,B) ;           \
    (B == NULL) ? GBI1_while (&Iter)   : GBI2_while (&Iter)    ;            \
    (B == NULL) ? GBI1_next  (&Iter)   : GBI2_next  (&Iter)

#define for_each_vector2m(A,B)                                              \
    GBI_iterator Iter ;                                                     \
    for (each_vector2m (Iter, A, B))

// iterate the union of three matrices A, B, and C
#define each_vector3(Iter,A,B,C)                                            \
    GBI3_init (&Iter,A,B,C) ; GBI3_while (&Iter)  ; GBI3_next (&Iter)

#define for_each_vector3(A,B,C)                                             \
    GBI_iterator Iter ;                                                     \
    for (each_vector3 (Iter, A, B, C))

// iterate over a matrix A and an implicit expanded scalar
#define each_vector2s(Iter,A)                                               \
    GBI2s_init (&Iter,A)    ; GBI2s_while (&Iter) ; GBI2s_next (&Iter)

#define for_each_vector2s(A)                                                \
    GBI_iterator Iter ;                                                     \
    for (each_vector2s (Iter, A))

// iterate over two matrices A and B, and an implicit expanded scalar
#define each_vector3s(Iter,A,B)                                             \
    GBI3s_init (&Iter,A,B)  ; GBI3s_while (&Iter) ; GBI3s_next (&Iter)

#define for_each_vector3s(A,B)                                              \
    GBI_iterator Iter ;                                                     \
    for (each_vector3s (Iter, A, B))

// iterate over a vector of a single matrix
#define each_entry(Iter,j,p,pend)                                           \
    int64_t GBI1_initj (Iter, j, p, pend) ; (p) < (pend) ; (p)++

#define for_each_entry(j,p,pend)                                            \
    int64_t GBI1_initj (Iter, j, p, pend) ;                                 \
    for ( ; (p) < (pend) ; (p)++)

//------------------------------------------------------------------------------
// GB_jstartup:  start the formation of a matrix
//------------------------------------------------------------------------------

// GB_jstartup is used with GB_jappend and GB_jwrapup to create the
// hyperlist and vector pointers of a new matrix, one at a time.

// GB_jstartup logs the start of C(:,0); it also acts as if it logs the end of
// the sentinal vector C(:,-1).

static inline void GB_jstartup
(
    GrB_Matrix C,           // matrix to start creating
    int64_t *jlast,         // last vector appended, set to -1
    int64_t *cnz,           // set to zero
    int64_t *cnz_last       // set to zero
)
{
    C->p [0] = 0 ;          // log the start of C(:,0)
    (*cnz) = 0 ;            //
    (*cnz_last) = 0 ;
    (*jlast) = -1 ;         // last sentinal vector is -1
    if (C->is_hyper)
    { 
        C->nvec = 0 ;       // clear all existing vectors from C
    }
    C->nvec_nonempty = 0 ;  // # of non-empty vectors will be counted
}

//------------------------------------------------------------------------------
// GB_jappend:  append a new vector to the end of a matrix
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2018, All Rights Reserved.
// http://suitesparse.com   See GraphBLAS/Doc/License.txt for license.

// Append a new vector to the end of a matrix C.

// If C->is_hyper is true, C is in hypersparse form with
// C->nvec <= C->plen <= C->vdim.  C->h has size C->plen.
// If C->is_hyper is false, C is in non-hypersparse form with
// C->nvec == C->plen == C->vdim.  C->h is NULL.
// In both cases, C->p has size C->plen+1.

// For both hypersparse and non-hypersparse, C->nvec_nonemty <= C->nvec
// is the number of vectors with at least one entry.

static inline GrB_Info GB_jappend
(
    GrB_Matrix C,           // matrix to append vector j to
    int64_t j,              // new vector to append
    int64_t *jlast,         // last vector appended, -1 if none
    int64_t cnz,            // nnz(C) after adding this vector j
    int64_t *cnz_last       // nnz(C) before adding this vector j
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    ASSERT (C != NULL) ;
    ASSERT (!C->p_shallow) ;
    ASSERT (!C->h_shallow) ;
    ASSERT (C->p != NULL) ;

    if (cnz <= (*cnz_last))
    { 
        // nothing to do
        return (GrB_SUCCESS) ;
    }

    // one more non-empty vector
    C->nvec_nonempty++ ;

    if (C->is_hyper)
    { 

        //----------------------------------------------------------------------
        // C is hypersparse; make sure space exists in the hyperlist
        //----------------------------------------------------------------------

        ASSERT (C->p [C->nvec] == (*cnz_last)) ;
        ASSERT (C->h != NULL) ;

        // check if space exists
        if (C->nvec == C->plen)
        { 
            // double the size of C->h and C->p
            GrB_Info info = GB_hyper_realloc (C, IMIN (C->vdim, 2*(C->plen+1)));
            if (info != GrB_SUCCESS)
            { 
                return (info) ;
            }
        }

        ASSERT (C->nvec >= 0) ;
        ASSERT (C->nvec < C->plen) ;
        ASSERT (C->plen <= C->vdim) ;
        ASSERT (C->p [C->nvec] == (*cnz_last)) ;

        C->h [C->nvec] = j ;            // add j to the hyperlist
        C->p [C->nvec+1] = cnz ;        // mark the end of C(:,j)
        C->nvec++ ;                     // one more vector in the hyperlist

    }
    else
    {

        //----------------------------------------------------------------------
        // C is non-hypersparse
        //----------------------------------------------------------------------

        int64_t *restrict Cp = C->p ;

        ASSERT (C->nvec == C->plen && C->plen == C->vdim) ;
        ASSERT (C->h == NULL) ;
        ASSERT (Cp [(*jlast)+1] == (*cnz_last)) ;

        // Even if C is non-hypersparse, the iteration that uses this function
        // may iterate over a hypersparse input matrix, so not every vector j
        // will be traversed.  So when j is seen, the end of vectors jlast+1 to
        // j must logged in Cp.

        for (int64_t jprior = (*jlast)+1 ; jprior < j ; jprior++)
        { 
            Cp [jprior+1] = (*cnz_last) ;   // mark the end of C(:,jprior)
        }
        Cp [j+1] = cnz ;                    // mark the end of C(:,j)
    }

    // record the last vector added to C
    (*cnz_last) = cnz ;
    (*jlast) = j ;

    return (GrB_SUCCESS) ;
}

//------------------------------------------------------------------------------
// GB_jwrapup:  finish contructing a new matrix
//------------------------------------------------------------------------------

// Log the end of any vectors in C that are not yet terminated.  Nothing
// happens if C is hypersparse (except for setting C->magic).

static inline void GB_jwrapup
(
    GrB_Matrix C,           // matrix to finish
    int64_t jlast,          // last vector appended, -1 if none
    int64_t cnz             // final nnz(C)
)
{

    if (!C->is_hyper)
    {

        //----------------------------------------------------------------------
        // C is non-hypersparse
        //----------------------------------------------------------------------

        // log the end of C(:,jlast+1) to C(:,n-1), in case the last vector
        // j=n-1 has not yet been seen, or has been seen but was empty.

        int64_t *restrict Cp = C->p ;
        int64_t j = C->vdim - 1 ;

        for (int64_t jprior = jlast+1 ; jprior <= j ; jprior++)
        { 
            Cp [jprior+1] = cnz ;           // mark the end of C(:,jprior)
        }
    }

    // C->p and C->h are now valid
    C->magic = MAGIC ;
}

#endif

