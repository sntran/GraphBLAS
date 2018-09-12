//------------------------------------------------------------------------------
// GB_stats: return memory usage and other statistics
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2018, All Rights Reserved.
// http://suitesparse.com   See GraphBLAS/Doc/License.txt for license.

//------------------------------------------------------------------------------

#include "GB.h"

void GB_stats
(
    GxB_Statistics *stats
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    ASSERT (stats != NULL) ;

    //--------------------------------------------------------------------------
    // get memory usage
    //--------------------------------------------------------------------------

    #pragma omp critical (GB_memory)
    {
        stats->nmalloc = GB_Global.nmalloc ;
        stats->inuse   = GB_Global.inuse ;
        stats->maxused = GB_Global.maxused ;
        GB_Global.maxused = GB_Global.inuse ;
    }

    //--------------------------------------------------------------------------
    // clear remainder of stats
    //--------------------------------------------------------------------------

    // these components are reserved for future use, so that new statistics can
    // be added without requiring a prior user application to be recompiled.

    for (int i = 0 ; i < 20 ; i++)
    { 
        stats->future [i] = 0 ;
        stats->xfuture [i] = 0 ;
    }
}

