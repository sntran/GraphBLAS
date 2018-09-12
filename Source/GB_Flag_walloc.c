//------------------------------------------------------------------------------
// GB_Flag_walloc: ensure Flag workspace is large enough
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2018, All Rights Reserved.
// http://suitesparse.com   See GraphBLAS/Doc/License.txt for license.

//------------------------------------------------------------------------------

#include "GB.h"

GrB_Info GB_Flag_walloc             // allocate Flag space
(
    int64_t Flag_required           // ensure Flag is at least this large
)
{ 

    int64_t currsize = GB_thread_local.Flag_size ;
    if (Flag_required > currsize)
    { 
        // free the existing space
        GB_FREE_MEMORY (GB_thread_local.Flag, currsize, sizeof (int8_t)) ;
        GB_thread_local.Flag_size = 0 ;

        // calloc the new space
        int64_t newsize = Flag_required + 1 ;
        GB_CALLOC_MEMORY (GB_thread_local.Flag, newsize, sizeof (int8_t)) ;
        if (GB_thread_local.Flag == NULL)
        { 
            // out of memory
            return (OUT_OF_MEMORY (GBYTES (newsize, sizeof (int8_t)))) ;
        }
        GB_thread_local.Flag_size = newsize ;
    }

    // this function can only be called when Flag [...] == 0
    // assertion for debugging only:
    ASSERT_FLAG_IS_CLEAR ;          // assert that Flag [...] == 0

    // success
    return (REPORT_SUCCESS) ;
}

