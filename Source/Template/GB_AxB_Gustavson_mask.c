//------------------------------------------------------------------------------
// GB_AxB_Gustavson_mask:  compute C<M>=A*B using the Gustavson method, with M
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2018, All Rights Reserved.
// http://suitesparse.com   See GraphBLAS/Doc/License.txt for license.

// This file is #include'd in GB_AxB_Gustavson.c, and Template/GB_AxB.c, the
// latter of which expands into Generated/GB_AxB__* for all built-in semirings.

// The pattern of C has not been computed, but NNZ(M) has given an upper bound
// on NNZ(C) so this method will not run out of memory.  This is Gustavson's
// method, extended to handle hypersparse matrices, arbitrary semirings, and a
// mask matrix M.

{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    ASSERT (NOT_ALIASED_3 (C, M, A, B)) ;
    ASSERT (C->vdim == B->vdim) ;
    ASSERT (C->vlen == A->vlen) ;
    ASSERT (A->vdim == B->vlen) ;

    ASSERT (C->vdim == M->vdim) ;
    ASSERT (C->vlen == M->vlen) ;

    //--------------------------------------------------------------------------
    // get workspace
    //--------------------------------------------------------------------------

    // get the Flag workspace (already allocated and cleared)
    int8_t *restrict Flag = GB_thread_local.Flag ;

    //--------------------------------------------------------------------------
    // get the mask
    //--------------------------------------------------------------------------

    const int64_t *restrict Mi = M->i ;
    const void    *restrict Mx = M->x ;
    GB_cast_function cast_M = GB_cast_factory (GB_BOOL_code, M->type->code) ;
    size_t msize = M->type->size ;

    //--------------------------------------------------------------------------
    // get A and B
    //--------------------------------------------------------------------------

    const int64_t *restrict Ap = A->p ;
    const int64_t *restrict Ai = A->i ;
    const int64_t *restrict Bi = B->i ;

    #ifdef HYPER
    const int64_t *restrict Ah = A->h ;
    int64_t anvec = A->nvec ;
    int64_t pleft, pright ;
    #endif

    //--------------------------------------------------------------------------
    // start the construction of C
    //--------------------------------------------------------------------------

    int64_t *restrict Ci = C->i ;

    int64_t jlast, cnz, cnz_last ;
    GB_jstartup (C, &jlast, &cnz, &cnz_last) ;

    //--------------------------------------------------------------------------
    // C<M>=A*B using the Gustavson method, pattern of C is a subset of M
    //--------------------------------------------------------------------------

    #ifdef HYPER
    for_each_vector2 (B, M)
    #else
    int64_t *restrict Bp = B->p ;
    int64_t *restrict Mp = M->p ;
    int64_t *restrict Cp = C->p ;
    int64_t n = C->vdim ;
    for (int64_t j = 0 ; j < n ; j++)
    #endif
    {

        //----------------------------------------------------------------------
        // get B(:,j) and M(:,j)
        //----------------------------------------------------------------------

        #ifdef HYPER
        int64_t GBI2_initj (Iter, j, pB_start, pB_end, pM_start, pM_end) ;
        #else
        int64_t pB_start = Bp [j] ;
        int64_t pB_end   = Bp [j+1] ;
        int64_t pM_start = Mp [j] ;
        int64_t pM_end   = Mp [j+1] ;
        #endif

        // C(:,j) is empty if either M(:,j) or B(:,j) are empty
        int64_t bjnz = pB_end - pB_start ;
        if (pM_start == pM_end || bjnz == 0)
        {
            #ifndef HYPER
            Cp [j+1] = cnz ;
            #endif
            continue ;
        }

        // M(:,j) has at least one entry; get the first and last index in M(:,j)
        int64_t im_first = Mi [pM_start] ;
        int64_t im_last  = Mi [pM_end-1] ;

        #ifdef HYPER
        // trim Ah on right
        if (A_is_hyper)
        { 
            pleft = 0 ;
            pright = anvec-1 ;
            if (bjnz > 2)
            {
                // trim Ah [0..pright] to remove any entries past last B(:,j)
                int64_t klast = Bi [pB_end-1] ;
                GB_bracket_right (klast, Ah, 0, &pright) ;
            }
        }
        #endif

        // M(:,j) is not yet scattered into Flag.  Flag is all zero
        bool marked = false ;

        //----------------------------------------------------------------------
        // C(:,j)<M(:,j)> = A * B(:,j), both values and pattern
        //----------------------------------------------------------------------

        for (int64_t pB = pB_start ; pB < pB_end ; pB++)
        {

            //------------------------------------------------------------------
            // get the pattern of B(k,j)
            //------------------------------------------------------------------

            int64_t k = Bi [pB] ;

            //------------------------------------------------------------------
            // get A(:,k)
            //------------------------------------------------------------------

            // find A(:,k), reusing pleft since Bi [...] is sorted
            int64_t pA_start, pA_end ;
            #ifdef HYPER
            GB_lookup (A_is_hyper, Ah, Ap, &pleft, pright, k,
                &pA_start, &pA_end) ;
            #else
            pA_start = Ap [k] ;
            pA_end   = Ap [k+1] ;
            #endif

            // skip if A(:,k) is empty
            if (pA_start == pA_end) continue ;

            // skip if the intersection of A(:,k) and M(:,j) is empty
            if (Ai [pA_end-1] < im_first || Ai [pA_start] > im_last) continue ;

            //------------------------------------------------------------------
            // scatter M(:,j) into Flag if not yet done
            //------------------------------------------------------------------

            if (!marked)
            {
                for (int64_t pM = pM_start ; pM < pM_end ; pM++)
                {
                    // mij = (bool) M (i,j)
                    bool mij ;
                    cast_M (&mij, Mx +(pM*msize), 0) ;
                    if (mij)
                    { 
                        // M(i,j) is true
                        Flag [Mi [pM]] = 1 ;
                    }
                }
                // M(:,j) has been scattered into Flag; must clear it when done
                marked = true ;

                // status of Flag [0..cvlen-1]:
                // Flag [i] = 0:  M(i,j) = 0, or entry not present in M(:,j)
                // Flag [i] = 1:  M(i,j) = 1,
            }

            //------------------------------------------------------------------
            // get the value of B(k,j)
            //------------------------------------------------------------------

            COPY_ARRAY_TO_SCALAR (bkj, Bx, pB, bsize) ;

            //------------------------------------------------------------------
            // w += (A(:,k) * B(k,j)) .* M(:,j)
            //------------------------------------------------------------------

            for (int64_t pA = pA_start ; pA < pA_end ; pA++)
            { 
                // w [i] += (A(i,k) * B(k,j)) .* M(i,j)
                int64_t i = Ai [pA] ;
                int8_t flag = Flag [i] ;
                if (flag == 0) continue ;
                // M(i,j) == 1 so do the work
                MULTADD_WITH_MASK ;
            }

            //------------------------------------------------------------------
            // status of Flag [0..cvlen-1] and w [0..cvlen-1]
            //------------------------------------------------------------------

            // Flag [i] = 0:  M(i,j)=0, or entry not present in M(:,j)
            // Flag [i] = 1:  M(i,j)=1, C(i,j) not present; w [i] uninitialized
            // Flag [i] = -1: M(i,j)=1, and C(i,j) is present; value is w [i]
        }

        //----------------------------------------------------------------------
        // check if C(:,j) is empty
        //----------------------------------------------------------------------

        // if M(:,j) has not been scattered into Flag, then C(:,j) must be
        // empty.  C(:,j) can still be empty if marked is false, but in that
        // case the Flag must still be cleared.

        #ifdef HYPER
        if (!marked) continue ;
        #endif

        //----------------------------------------------------------------------
        // gather C(:,j), both values and pattern, from pattern of M(:,j) and w
        //----------------------------------------------------------------------

        if (marked)
        {
            for (int64_t pM = pM_start ; pM < pM_end ; pM++)
            {
                int64_t i = Mi [pM] ;
                if (Flag [i] < 0)
                { 
                    // C(i,j) is a live entry, gather its row and value
                    // Cx [cnz] = w [i] ;
                    COPY_ARRAY_TO_ARRAY (Cx, cnz, w, i, zsize) ;
                    Ci [cnz++] = i ;
                }
                // clear the Flag
                Flag [i] = 0 ;
            }
        }

        //----------------------------------------------------------------------
        // log the end of C(:,j)
        //----------------------------------------------------------------------

        #ifdef HYPER
        // cannot fail since C->plen is the upper bound: number of non-empty
        // columns of B
        info = GB_jappend (C, j, &jlast, cnz, &cnz_last) ;
        ASSERT (info == GrB_SUCCESS) ;
        #else
        Cp [j+1] = cnz ;
        if (cnz > cnz_last) C->nvec_nonempty++ ;
        cnz_last = cnz ;
        #endif
    }

    //--------------------------------------------------------------------------
    // finalize C
    //--------------------------------------------------------------------------

    #ifdef HYPER
    GB_jwrapup (C, jlast, cnz) ;
    ASSERT (info == GrB_SUCCESS) ;
    #else
    C->magic = MAGIC ;
    #endif
}

