//------------------------------------------------------------------------------
// GxB_stats: return memory usage and other statistics
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2018, All Rights Reserved.
// http://suitesparse.com   See GraphBLAS/Doc/License.txt for license.

//------------------------------------------------------------------------------

#include "GB.h"

GrB_Info GxB_stats
(
    GxB_Statistics *stats
)
{ 

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    WHERE ("GxB_stats (&stats) ;") ;
    RETURN_IF_NULL (stats) ;

    //--------------------------------------------------------------------------
    // get statistics
    //--------------------------------------------------------------------------

    GB_stats (stats) ;
    return (REPORT_SUCCESS) ;
}

