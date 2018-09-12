//------------------------------------------------------------------------------
// GB_free: free a matrix
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2018, All Rights Reserved.
// http://suitesparse.com   See GraphBLAS/Doc/License.txt for license.

//------------------------------------------------------------------------------

// free all the content of a matrix.  After GB_free (&A), A is set to NULL.

#include "GB.h"

void GB_free                    // free a matrix
(
    GrB_Matrix *matrix          // handle of matrix to free
)
{

    if (matrix != NULL)
    {
        GrB_Matrix A = *matrix ;
        if (A != NULL && (A->magic == MAGIC || A->magic == MAGIC2))
        { 
            // free all content of A
            GB_phix_free (A) ;
            // free the header of A itself
            A->magic = FREED ;      // to help detect dangling pointers
            GB_FREE_MEMORY (*matrix, 1, sizeof (struct GB_Matrix_opaque)) ;
        }
        (*matrix) = NULL ;
    }
}

