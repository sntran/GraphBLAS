//------------------------------------------------------------------------------
// GB_AxB_heap_meta: compute C<M>=A*B or C=A*B using a heap-based method
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2018, All Rights Reserved.
// http://suitesparse.com   See GraphBLAS/Doc/License.txt for license.

//------------------------------------------------------------------------------

{
    if (M != NULL)
    { 
        // C<M> = A*B via a heap
        #define MASK
        #include "GB_AxB_heap_mask.c"
        #undef MASK
    }
    else
    { 
        // C = A*B via the heap
        #include "GB_AxB_heap_mask.c"
    }
}

