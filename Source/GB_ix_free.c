//------------------------------------------------------------------------------
// GB_ix_free: free A->i, A->x, pending tuples, zombies; A->p, A->h unchanged
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2018, All Rights Reserved.
// http://suitesparse.com   See GraphBLAS/Doc/License.txt for license.

//------------------------------------------------------------------------------

// Since A->p and A->h are unchanged, the matrix is still valid (unless it was
// invalid on input).  NNZ(A) would report zero, and so would GrB_Matrix_nvals.

#include "GB.h"

void GB_ix_free                 // free A->i and A->x of a matrix
(
    GrB_Matrix A                // matrix with content to free
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    if (A == NULL)
    { 
        return ;
    }

    //--------------------------------------------------------------------------
    // free all but A->p and A->h
    //--------------------------------------------------------------------------

    // zombies and pending tuples are about to be deleted
    ASSERT (PENDING_OK (A)) ;
    ASSERT (ZOMBIES_OK (A)) ;

    // free A->i unless it is shallow
    if (!A->i_shallow) GB_FREE_MEMORY (A->i, A->nzmax, sizeof (int64_t)) ;
    A->i = NULL ;
    A->i_shallow = false ;

    // free A->x unless it is shallow
    if (!A->x_shallow) GB_FREE_MEMORY (A->x, A->nzmax, A->type->size) ;
    A->x = NULL ;
    A->x_shallow = false ;

    A->nzmax = 0 ;

    // no zombies remain
    A->nzombies = 0 ;

    // free pending tuples
    GB_pending_free (A) ;

    // remove from the queue, if present
    GB_queue_remove (A) ;
}

