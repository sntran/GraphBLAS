//------------------------------------------------------------------------------
// GB_global_option_set: set one or more global options
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2018, All Rights Reserved.
// http://suitesparse.com   See GraphBLAS/Doc/License.txt for license.

//------------------------------------------------------------------------------

#include "GB.h"

void GB_global_option_set       // set one or more global options
(
    bool set_hyper_ratio,       // if true, set the hyper_ratio
    double hyper_ratio,
    bool set_format,            // if true, set the format
    GxB_Format_Value format
)
{

    // split into multiple critical sections since it is likely the writes
    // are already atomic, so no synchronization needs to be done.

    if (set_hyper_ratio)
    { 
        #pragma omp critical (GB_options)
        {
            GB_Global.hyper_ratio = hyper_ratio ;
        }
    }

    if (set_format)
    { 
        bool is_csc = (format != GxB_BY_ROW) ;
        #pragma omp critical (GB_options)
        {
            GB_Global.is_csc = is_csc ;
        }
    }
}

