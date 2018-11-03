//------------------------------------------------------------------------------
// GrB_finalize: finalize GraphBLAS
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2018, All Rights Reserved.
// http://suitesparse.com   See GraphBLAS/Doc/License.txt for license.

//------------------------------------------------------------------------------

// GrB_finalize must be called as the last GraphBLAS function.

// In this version of SuiteSparse:GraphBLAS, GrB_finalize frees the workspace
// held internally in thread-local storage.  It can be called at any time and
// can be followed by GraphBLAS function.

// The error condition is not modified with "return (GB_REPORT_SUCCESS)" so that
// GrB_error can be called to report any prior error.

#include "GB.h"

GrB_Info GrB_finalize ( )
{ 

    GB_wfree ( ) ;              // free all thread-local workspace
    return (GrB_SUCCESS) ;      // method always succeeds
}

