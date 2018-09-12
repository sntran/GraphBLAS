//------------------------------------------------------------------------------
// GB_Descriptor_check: check and print a Descriptor
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2018, All Rights Reserved.
// http://suitesparse.com   See GraphBLAS/Doc/License.txt for license.

//------------------------------------------------------------------------------

#include "GB.h"

//------------------------------------------------------------------------------
// dcheck: check a single descriptor field
//------------------------------------------------------------------------------

static void dcheck
(
    GrB_Info *info,
    bool spec,
    const char *field,
    const GrB_Desc_Value v,
    const GrB_Desc_Value nondefault,
    int pr
)
{

    bool ok = true ;

    if (pr > 0) printf ("D.%s = ", field) ;
    switch (v)
    {
        case GxB_DEFAULT       : if (pr > 0) printf ("default   ") ; break ;
        case GrB_SCMP          : if (pr > 0) printf ("scmp      ") ; break ;
        case GrB_TRAN          : if (pr > 0) printf ("tran      ") ; break ;
        case GrB_REPLACE       : if (pr > 0) printf ("replace   ") ; break ;
        case GxB_AxB_GUSTAVSON : if (pr > 0) printf ("Gustavson ") ; break ;
        case GxB_AxB_HEAP      : if (pr > 0) printf ("heap      ") ; break ;
        case GxB_AxB_DOT       : if (pr > 0) printf ("dot       ") ; break ;
        default                : if (pr > 0) printf ("unknown   ") ;
            *info = GrB_INVALID_OBJECT ;
            ok = false ;
            break ;
    }

    if (ok)
    { 
        if (spec)
        { 
            // descriptor field can be set to the default,
            // or one non-default value
            if (! (v == GxB_DEFAULT || v == nondefault))
            { 
                ok = false ;
            }
        }
        else
        { 
            // GxB_AxB_METHOD:
            if (! (v == GxB_DEFAULT || v == GxB_AxB_GUSTAVSON
                || v == GxB_AxB_HEAP || v == GxB_AxB_DOT))
            { 
                ok = false ;
            }
        }
    }

    if (!ok)
    { 
        if (pr > 0) printf (" (invalid value for this field)") ;
        *info = GrB_INVALID_OBJECT ;
    }

    if (pr > 0) printf ("\n") ;
}

//------------------------------------------------------------------------------
// GB_Descriptor_check
//------------------------------------------------------------------------------

GrB_Info GB_Descriptor_check    // check a GraphBLAS descriptor
(
    const GrB_Descriptor D,     // GraphBLAS descriptor to print and check
    const char *name,           // name of the descriptor, optional
    int pr                      // 0: print nothing, 1: print header and
                                // errors, 2: print brief, 3: print all
)
{ 

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    if (pr > 0) printf ("\nGraphBLAS Descriptor: %s ", NAME) ;

    if (D == NULL)
    { 
        // GrB_error status not modified since this may be an optional argument
        if (pr > 0) printf ("NULL\n") ;
        return (GrB_NULL_POINTER) ;
    }

    //--------------------------------------------------------------------------
    // check object
    //--------------------------------------------------------------------------

    CHECK_MAGIC (D, "Descriptor") ;

    if (pr > 0) printf ("\n") ;

    GrB_Info info = GrB_SUCCESS ;
    dcheck (&info, true,  "output    ", D->out,  GrB_REPLACE, pr) ;
    dcheck (&info, true,  "mask      ", D->mask, GrB_SCMP,    pr) ;
    dcheck (&info, true,  "input0    ", D->in0,  GrB_TRAN,    pr) ;
    dcheck (&info, true,  "input1    ", D->in1,  GrB_TRAN,    pr) ;
    dcheck (&info, false, "AxB_method", D->axb,  0,           pr) ;

    if (info != GrB_SUCCESS)
    { 
        if (pr > 0) printf ("Descriptor field set to an invalid value\n") ;
        return (ERROR (GrB_INVALID_OBJECT, (LOG,
            "Descriptor field set to an invalid value: [%s]", NAME))) ;
    }

    return (GrB_SUCCESS) ; // not REPORT_SUCCESS; may mask error in caller
}

