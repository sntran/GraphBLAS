//------------------------------------------------------------------------------
// GB_global_option_get: get one or more global options
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2018, All Rights Reserved.
// http://suitesparse.com   See GraphBLAS/Doc/License.txt for license.

//------------------------------------------------------------------------------

#include "GB.h"

void GB_global_option_get       // get one or more global options
(
    double *hyper_ratio,
    GxB_Format_Value *format,
    bool *is_csc
)
{

    // split into multiple critical sections since it is likely the reads
    // are already atomic, so no synchronization needs to be done.

    if (hyper_ratio != NULL)
    { 
        double h ;
        #pragma omp critical (GB_options)
        {
            h = GB_Global.hyper_ratio ;
        }
        (*hyper_ratio) = h ;
    }

    if (format != NULL || is_csc != NULL)
    { 
        bool c ;
        #pragma omp critical (GB_options)
        {
            c = GB_Global.is_csc ;
        }
        if (is_csc != NULL)
        { 
            (*is_csc) = c ;
        }
        if (format != NULL)
        { 
            (*format) = (c) ? GxB_BY_COL : GxB_BY_ROW ;
        }
    }
}

