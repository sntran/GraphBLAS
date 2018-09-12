//------------------------------------------------------------------------------
// GrB_Vector_setElement: set an entry in a vector, w (row) = x
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2018, All Rights Reserved.
// http://suitesparse.com   See GraphBLAS/Doc/License.txt for license.

//------------------------------------------------------------------------------

// Set a single scalar, w(row) = x, typecasting from the type of x to
// the type of w as needed.

#include "GB.h"

#define SET(type,T,ampersand)                                               \
GrB_Info GrB_Vector_setElement_ ## T    /* w(row) = x    */                 \
(                                                                           \
    GrB_Vector w,                       /* vector to modify           */    \
    const type x,                       /* scalar to assign to w(row) */    \
    GrB_Index row                       /* row index                  */    \
)                                                                           \
{                                                                           \
    WHERE ("GrB_Vector_setElement_" GB_STR(T) " (w, x, row)") ;             \
    RETURN_IF_NULL_OR_FAULTY (w) ;                                          \
    ASSERT (VECTOR_OK (w)) ;                                                \
    return (GB_setElement ((GrB_Matrix) w, ampersand x, row, 0,             \
        GB_ ## T ## _code)) ;                                               \
}

SET (bool     , BOOL   , &) ;
SET (int8_t   , INT8   , &) ;
SET (uint8_t  , UINT8  , &) ;
SET (int16_t  , INT16  , &) ;
SET (uint16_t , UINT16 , &) ;
SET (int32_t  , INT32  , &) ;
SET (uint32_t , UINT32 , &) ;
SET (int64_t  , INT64  , &) ;
SET (uint64_t , UINT64 , &) ;
SET (float    , FP32   , &) ;
SET (double   , FP64   , &) ;
SET (void *   , UDT    ,  ) ;

#undef SET

