//------------------------------------------------------------------------------
// GxB_Global_Option_set: set a global default option for all future matrices
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2018, All Rights Reserved.
// http://suitesparse.com   See GraphBLAS/Doc/License.txt for license.

//------------------------------------------------------------------------------

#include "GB.h"

GrB_Info GxB_Global_Option_set      // set a global default option
(
    const GxB_Option_Field field,   // option to change
    ...                             // value to change it to
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    GB_WHERE ("GxB_Global_Option_set (field, value)") ;

    //--------------------------------------------------------------------------
    // set the global option
    //--------------------------------------------------------------------------

    va_list ap ;
    double hyper_ratio ;
    GxB_Format_Value format ;

    switch (field)
    {

        case GxB_HYPER : 

            va_start (ap, field) ;
            hyper_ratio = va_arg (ap, double) ;
            va_end (ap) ;
            GB_global_option_set (true, hyper_ratio, false, 0) ;
            break ;

        case GxB_FORMAT : 

            va_start (ap, field) ;
            format = va_arg (ap, GxB_Format_Value) ;
            va_end (ap) ;
            GB_global_option_set (false, 0, true, format) ;
            break ;

        default : 

            return (GB_ERROR (GrB_INVALID_VALUE, (GB_LOG,
                    "invalid option field [%d], must be one of:\n"
                    "GxB_HYPER [%d] or GxB_FORMAT [%d]",
                    field, GxB_HYPER, GxB_FORMAT))) ;

    }

    return (GB_REPORT_SUCCESS) ;
}

