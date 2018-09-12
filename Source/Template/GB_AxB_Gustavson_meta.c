//------------------------------------------------------------------------------
// GB_AxB_Gustavson_meta: C=A*B and C<M>=A*B
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2018, All Rights Reserved.
// http://suitesparse.com   See GraphBLAS/Doc/License.txt for license.

//------------------------------------------------------------------------------

{
    bool A_is_hyper = IS_HYPER (A) ;
    if (A_is_hyper || IS_HYPER (B) || IS_HYPER (C) || IS_HYPER (M))
    {
        #define HYPER
        if (M != NULL)
        { 
            // C<M> = A*B where M is pattern of C
            #include "GB_AxB_Gustavson_mask.c"
        }
        else
        { 
            // C = A*B with pattern of C as defined on input
            #include "GB_AxB_Gustavson_nomask.c"
        }
        #undef HYPER
    }
    else
    {
        if (M != NULL)
        { 
            // C<M> = A*B where M is pattern of C
            #include "GB_AxB_Gustavson_mask.c"
        }
        else
        { 
            // C = A*B with pattern of C as defined on input
            #include "GB_AxB_Gustavson_nomask.c"
        }
    }
}
