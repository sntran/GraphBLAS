//------------------------------------------------------------------------------
// GrB_Matrix_setElement: set an entry in a matrix, C(row,col) = x
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2018, All Rights Reserved.
// http://suitesparse.com   See GraphBLAS/Doc/License.txt for license.

//------------------------------------------------------------------------------

// Set a single entry in a matrix, C(row,col) = x in MATLAB notation,
// typecasting from the type of x to the type of C, as needed.

#include "GB.h"

#define SET(type,T,ampersand)                                               \
GrB_Info GrB_Matrix_setElement_ ## T    /* C (row,col) = x */               \
(                                                                           \
    GrB_Matrix C,                       /* matrix to modify               */\
    const type x,                       /* scalar to assign to C(row,col) */\
    GrB_Index row,                      /* row index                      */\
    GrB_Index col                       /* column index                   */\
)                                                                           \
{                                                                           \
    WHERE ("GrB_Matrix_setElement_" GB_STR(T) " (C, row, col, x)") ;        \
    RETURN_IF_NULL_OR_FAULTY (C) ;                                          \
    return (GB_setElement (C, ampersand x, row, col, GB_ ## T ## _code)) ;  \
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

