//------------------------------------------------------------------------------
// GxB_Global_Option_get: get a global default option for all future matrices
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2018, All Rights Reserved.
// http://suitesparse.com   See GraphBLAS/Doc/License.txt for license.

//------------------------------------------------------------------------------

#include "GB.h"

GrB_Info GxB_Global_Option_get      // gets the current global option
(
    const GxB_Option_Field field,   // option to query
    ...                             // return value of the global option
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    WHERE ("GxB_Global_Option_get (field, &value)") ;

    //--------------------------------------------------------------------------
    // get the option
    //--------------------------------------------------------------------------

    va_list ap ;
    double *hyper_ratio, hyper ;
    GxB_Format_Value *format, fmt ;

    GB_global_option_get (&hyper, &fmt, NULL) ;

    switch (field)
    {

        case GxB_HYPER : 

            va_start (ap, field) ;
            hyper_ratio = va_arg (ap, double *) ;
            va_end (ap) ;

            RETURN_IF_NULL (hyper_ratio) ;
            (*hyper_ratio) = hyper ;
            break ;

        case GxB_FORMAT : 

            va_start (ap, field) ;
            format = va_arg (ap, GxB_Format_Value *) ;
            va_end (ap) ;

            RETURN_IF_NULL (format) ;
            (*format) = fmt ;
            break ;

        default : 

            return (ERROR (GrB_INVALID_VALUE, (LOG,
                    "invalid option field [%d], must be one of:\n"
                    "GxB_HYPER [%d] or GxB_FORMAT [%d]",
                    field, GxB_HYPER, GxB_FORMAT))) ;

    }

    return (REPORT_SUCCESS) ;
}

