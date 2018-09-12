//------------------------------------------------------------------------------
// GB_SelectOp_new: create a new select operator
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2018, All Rights Reserved.
// http://suitesparse.com   See GraphBLAS/Doc/License.txt for license.

//------------------------------------------------------------------------------

// The select function signature must be:

//      bool f (GrB_Index i, GrB_Index j, GrB_Index nrows, GrB_Index ncols,
//              const void *x, const void *k) ;

// This function is not directly user-callable.  Use GxB_SelectOp_new instead.

#include "GB.h"

GrB_Info GB_SelectOp_new        // create a new user-defined select operator
(
    GxB_SelectOp *selectop,     // handle for the new select operator
    void *function,             // pointer to the select function
    const GrB_Type xtype,       // type of input x
    const char *name            // name of the function
)
{ 

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    WHERE ("GxB_SelectOp_new (selectop, function, xtype)") ;
    RETURN_IF_NULL (selectop) ;
    (*selectop) = NULL ;
    RETURN_IF_NULL (function) ;
    RETURN_IF_FAULTY (xtype) ;   // xtype may be NULL

    //--------------------------------------------------------------------------
    // create the select op
    //--------------------------------------------------------------------------

    // allocate the select operator
    GB_CALLOC_MEMORY (*selectop, 1, sizeof (struct GB_SelectOp_opaque)) ;
    if (*selectop == NULL)
    { 
        return (NO_MEMORY) ;
    }

    // initialize the select operator
    GxB_SelectOp op = *selectop ;
    op->magic = MAGIC ;
    op->xtype = xtype ;
    op->function = function ;
    strncpy (op->name, name, GB_LEN-1) ;
    op->opcode = GB_USER_SELECT_opcode ;    // generic opcode for all user ops
    ASSERT_OK (GB_check (op, "new user-defined select op", D0)) ;
    return (REPORT_SUCCESS) ;
}

