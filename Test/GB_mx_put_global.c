//------------------------------------------------------------------------------
// GB_mx_put_global: put the GraphBLAS status in MATLAB workspace
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2018, All Rights Reserved.
// http://suitesparse.com   See GraphBLAS/Doc/License.txt for license.

//------------------------------------------------------------------------------

#include "GB_mex.h"

void GB_mx_put_global
(
    bool cover
)
{

    //--------------------------------------------------------------------------
    // check nmalloc
    //--------------------------------------------------------------------------

    Complex_finalize ( ) ;

    //--------------------------------------------------------------------------
    // return the time to MATLAB, if it was computed
    //--------------------------------------------------------------------------

    GB_mx_put_time ( ) ;

    //--------------------------------------------------------------------------
    // log the statement coverage
    //--------------------------------------------------------------------------

    #ifdef GBCOVER
    if (cover) GB_cover_put ( ) ;
    #endif

    //--------------------------------------------------------------------------
    // finalize GraphBLAS
    //--------------------------------------------------------------------------

    GB_wfree ( ) ;

    GxB_Statistics stats ;
    GB_stats (&stats) ;
    if (stats.nmalloc != 0)
    {
        printf ("GraphBLAS nmalloc "GBd"! inuse "GBd" maxused "GBd"\n",
            stats.nmalloc, stats.inuse, stats.maxused) ;
        mexErrMsgTxt ("memory leak!") ;
    }
}

