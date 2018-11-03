//------------------------------------------------------------------------------
// GB_Work_walloc: ensure Work workspace is large enough
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2018, All Rights Reserved.
// http://suitesparse.com   See GraphBLAS/Doc/License.txt for license.

//------------------------------------------------------------------------------

#include "GB.h"

GrB_Info GB_Work_walloc             // allocate Work space
(
    size_t Work_nitems_required,    // ensure Work is at least this large,
    size_t Work_itemsize            // # items times size of each item
)
{

    size_t Work_required ;

    // Work_required = Work_nitems_required * Work_itemsize
    if (!GB_size_t_multiply (&Work_required,
        Work_nitems_required, Work_itemsize))
    { 
        // size_t overflow
        double memory = GBYTES (Work_nitems_required, Work_itemsize) ;
        return (GB_OUT_OF_MEMORY (memory)) ;
    }

    if (Work_required > GB_thread_local.Work_size)
    {
        // free the existing space
        GB_Work_wfree ( ) ;

        // malloc the new space
        size_t newsize = Work_required + 1 ;
        GB_MALLOC_MEMORY (GB_thread_local.Work, newsize, sizeof (char)) ;
        if (GB_thread_local.Work == NULL)
        { 
            // out of memory
            return (GB_OUT_OF_MEMORY (GBYTES (newsize, sizeof (char)))) ;
        }
        GB_thread_local.Work_size = newsize ;
    }

    return (GB_REPORT_SUCCESS) ;
}

