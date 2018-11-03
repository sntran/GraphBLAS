//------------------------------------------------------------------------------
// GxB_stats: return memory usage and other statistics
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2018, All Rights Reserved.
// http://suitesparse.com   See GraphBLAS/Doc/License.txt for license.

//------------------------------------------------------------------------------

// This function has never appeared in the user guide and will be deprecated.
// Its functionality will be replaced by GxB_get.  The function will be deleted
// in SuiteSparse:GraphBLAS 3.0.

#include "GB.h"

GrB_Info GxB_stats
(
    GxB_Statistics *stats
)
{ 

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    GB_WHERE ("GxB_stats (&stats) ;") ;
    GB_RETURN_IF_NULL (stats) ;

    //--------------------------------------------------------------------------
    // get statistics
    //--------------------------------------------------------------------------

    GB_stats (stats) ;
    return (GB_REPORT_SUCCESS) ;
}

