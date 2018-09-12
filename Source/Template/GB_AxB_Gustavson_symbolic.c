//------------------------------------------------------------------------------
// GB_AxB_Gustavson_symbolic: C=A*B symbolic analysis
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2018, All Rights Reserved.
// http://suitesparse.com   See GraphBLAS/Doc/License.txt for license.

//------------------------------------------------------------------------------

{

    //--------------------------------------------------------------------------
    // get A and B
    //--------------------------------------------------------------------------

    const int64_t *Bp = B->p ;
    const int64_t *Bi = B->i ;
    const int64_t *Ap = A->p ;
    const int64_t *Ai = A->i ;

    #ifdef HYPER
    const int64_t *Ah = A->h ;
    int64_t anvec = A->nvec ;
    #endif

    //--------------------------------------------------------------------------
    // start the construction of the pattern of C
    //--------------------------------------------------------------------------

    int64_t *restrict Ci = C->i ;
    int64_t *restrict Cp = C->p ;

    int64_t jlast, cnz, cnz_last ;
    GB_jstartup (C, &jlast, &cnz, &cnz_last) ;

    //--------------------------------------------------------------------------
    // symbolic pattern of C = A*B
    //--------------------------------------------------------------------------

    #ifdef HYPER
    for_each_vector (B)
    #else
    int64_t n = C->vdim ;
    for (int64_t j = 0 ; j < n ; j++)
    #endif
    {

        #ifdef HYPER
        int64_t GBI1_initj (Iter, j, pB_start, pB_end) ;
        #else
        int64_t pB_start = Bp [j] ;
        int64_t pB_end   = Bp [j+1] ;
        #endif

        //----------------------------------------------------------------------
        // reallocate C if necessary
        //----------------------------------------------------------------------

        // Note that cvlen is an upper bound on nnz (C (:,j)), but it can
        // be a very loose bound if C is hypersparse.
        int64_t cmax = cnz + cvlen ;
        if (cmax > C->nzmax)
        { 
            OK (GB_ix_realloc (C, 4*(C->nzmax + cvlen), false)) ;
            Ci = C->i ;
        }

        //----------------------------------------------------------------------
        // C(:,j) = set union of all A(:,k) for each nonzero B(k,j) ;
        //----------------------------------------------------------------------

        int64_t bjnz = pB_end - pB_start ;

        #ifdef HYPER
        int64_t pleft = 0 ;
        int64_t pright = anvec-1 ;
        #endif

        if (bjnz == 0)
        { 

            //------------------------------------------------------------------
            // B (:,j) is empty; nothing to do
            //------------------------------------------------------------------

            #ifdef HYPER
            continue ;
            #endif

        }
        else if (bjnz == 1)
        {

            //------------------------------------------------------------------
            // C (:,j) = A (:,k) for a single nonzero B(k,j)
            //------------------------------------------------------------------

            // C(:,j) = A(:,k)
            int64_t k = Bi [pB_start] ;

            // find A(:,k)
            int64_t pA_start, pA_end ;
            #ifdef HYPER
            GB_lookup (A_is_hyper, Ah, Ap, &pleft, pright, k,
                &pA_start, &pA_end) ;
            #else
            pA_start = Ap [k] ;
            pA_end   = Ap [k+1] ;
            #endif

            for (int64_t pA = pA_start ; pA < pA_end ; pA++)
            { 
                int64_t i = Ai [pA] ;
                // C(i,j) is nonzero
                Ci [cnz++] = i ;
            }

        }
        else if (bjnz == 2)
        {

            //------------------------------------------------------------------
            // 2-way merge of A (:,k1) and A (:,k2)
            //------------------------------------------------------------------

            int64_t k1 = Bi [pB_start] ;
            int64_t k2 = Bi [pB_start+1] ;
            ASSERT (k1 < k2) ;

            int64_t p1, p1_end, p2, p2_end ;

            // find A(:,k1) and A(:,k2)
            #ifdef HYPER
            GB_lookup (A_is_hyper, Ah, Ap, &pleft, pright, k1,
                &p1, &p1_end) ;
            // Use pleft of k1 to trim the search for k2 since k1 < k2
            GB_lookup (A_is_hyper, Ah, Ap, &pleft, pright, k2,
                &p2, &p2_end) ;
            #else
            p1     = Ap [k1] ;
            p1_end = Ap [k1+1] ;
            p2     = Ap [k2] ;
            p2_end = Ap [k2+1] ;
            #endif

            while (p1 < p1_end || p2 < p2_end)
            {
                int64_t i1 = (p1 < p1_end) ? Ai [p1] : cvlen ;
                int64_t i2 = (p2 < p2_end) ? Ai [p2] : cvlen ;
                int64_t i ;
                if (i1 < i2)
                { 
                    i = i1 ;
                    p1++ ;
                }
                else if (i1 > i2)
                { 
                    i = i2 ;
                    p2++ ;
                }
                else // i1 == i2
                { 
                    i = i1 ;
                    p1++ ;
                    p2++ ;
                }
                // C(i,j) is nonzero
                Ci [cnz++] = i ;
            }

        }
        else
        {

            //------------------------------------------------------------------
            // general case, nnz (B (:,j)) > 2
            //------------------------------------------------------------------

            // flag++
            int64_t flag = GB_Mark_reset (1, 0) ;

            #ifdef HYPER
            // trim on right
            if (A_is_hyper)
            { 
                // trim Ah [0..pright] to remove any entries past
                // the last B(:,j)
                GB_bracket_right (Bi [pB_end-1], Ah, 0, &pright) ;
            }
            #endif

            for (int64_t pB = pB_start ; pB < pB_end ; pB++)
            {
                // if C(:,j) now completely full, no need to continue
                if (cnz == cmax) break ;

                // symbolic saxpy C(:,j) += A(:,k)*B(k,j)
                int64_t k = Bi [pB] ;

                // find A(:,k), reusing pleft since Bi [...] is sorted
                int64_t pA_start, pA_end ;
                #ifdef HYPER
                GB_lookup (A_is_hyper, Ah, Ap, &pleft, pright, k,
                    &pA_start, &pA_end) ;
                #else
                pA_start = Ap [k] ;
                pA_end   = Ap [k+1] ;
                #endif

                for (int64_t pA = pA_start ; pA < pA_end ; pA++)
                {
                    int64_t i = Ai [pA] ;
                    // C(i,j) is nonzero
                    if (Mark [i] < flag)
                    { 
                        // C(i,j) is nonzero, and this is the 1st time row
                        // i has been added to the pattern in C(:,j).  Mark
                        // it so row i is not added again.
                        Mark [i] = flag ;
                        // add to the column pattern of A*B
                        Ci [cnz++] = i ;
                    }
                }
            }

            // sort the pattern of C(:,j)
            int64_t len = cnz - cnz_last ;
            if (len == cvlen)
            {
                // no need to sort C(:,j) if dense; just recreate it
                for (int64_t pC = cnz_last, i = 0 ; pC < cnz ; pC++, i++)
                { 
                    Ci [pC] = i ;
                }
            }
            else
            { 
                // sort the nonzero indices in C(:,j)
                GB_qsort_1 (Ci + cnz_last, len) ;
            }
        }

        //----------------------------------------------------------------------
        // log the end of vector C(:,j)
        //----------------------------------------------------------------------

        #ifdef HYPER
        // this cannot fail since C->plen is the upper bound: the number
        // of non-empty columns of B.
        info = (GB_jappend (C, j, &jlast, cnz, &cnz_last)) ;
        ASSERT (info == GrB_SUCCESS) ;
        // if it could fail:
        // OK (info) ;              // check result and return on error
        #else
        Cp [j+1] = cnz ;
        if (cnz > cnz_last) C->nvec_nonempty++ ;
        cnz_last = cnz ;
        #endif

        // it also cannot run out of space here, but can do so above
        ASSERT (cnz <= C->nzmax) ;
    }

    //--------------------------------------------------------------------------
    // finalize C and clear the Mark
    //--------------------------------------------------------------------------

    #ifdef HYPER
    GB_jwrapup (C, jlast, cnz) ;
    #else
    C->magic = MAGIC ;
    #endif

    // clear the Mark array
    GB_Mark_reset (1, 0) ;

    //--------------------------------------------------------------------------
    // reduce the size of C->i to hold just the required space
    //--------------------------------------------------------------------------

    info = GB_ix_realloc (C, cnz, false) ;
    ASSERT (info == GrB_SUCCESS) ;
    ASSERT_OK (GB_check (C, "C symbolic Gustavson C=A*B", D0)) ;
}
