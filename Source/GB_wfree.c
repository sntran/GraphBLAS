//------------------------------------------------------------------------------
// GB_wfree: free all workspace
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2018, All Rights Reserved.
// http://suitesparse.com   See GraphBLAS/Doc/License.txt for license.

//------------------------------------------------------------------------------

// Frees all workspace held internally in thread-local storage.  It can be
// called at any time and can be followed by GraphBLAS function.

#include "GB.h"

void GB_wfree ( )
{ 

    // free all thread-local workspace
    GB_Mark_wfree ( ) ;
    GB_Work_wfree ( ) ;
    GB_Flag_wfree ( ) ;
}

