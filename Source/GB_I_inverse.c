//------------------------------------------------------------------------------
// GB_I_inverse: invert an index list
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2018, All Rights Reserved.
// http://suitesparse.com   See GraphBLAS/Doc/License.txt for license.

//------------------------------------------------------------------------------

// I is a large list relative to the vector length, avlen, and it is not
// contiguous.  Scatter I into the I inverse buckets (Mark and Inext) for quick
// lookup.

#include "GB.h"

GrB_Info GB_I_inverse
(
    const GrB_Index *I,         // list of indices, duplicates OK
    int64_t nI,                 // length of I
    int64_t avlen,              // length of the vectors of A
    bool need_Iwork1,           // true if Iwork1 of size nI needed for sorting
    // outputs:
    int64_t **p_Mark,           // head pointers for buckets, size avlen
    int64_t **p_Inext,          // next pointers for buckets, size nI
    int64_t **p_Iwork1,         // workspace of size nI, if needed
    int64_t *p_nduplicates,     // number of duplicate entires in I
    int64_t *p_flag             // Mark [0:avlen-1] < flag
)
{

    //--------------------------------------------------------------------------
    // get inputs
    //--------------------------------------------------------------------------

    int64_t *Mark = NULL ;
    int64_t *Inext = NULL ;
    int64_t *Iwork1 = NULL ;
    int64_t nduplicates = 0 ;
    int64_t flag = 0 ;

    *p_Mark = NULL ;
    *p_Inext = NULL ;
    *p_Iwork1 = NULL ;
    *p_nduplicates = 0 ;
    *p_flag = 0 ;

    //--------------------------------------------------------------------------
    // ensure Work is large enough for the scattered form of I
    //--------------------------------------------------------------------------

    int64_t iworksize = nI ;            // Inext [nI]

    if (need_Iwork1)
    { 
        iworksize += nI ;               // Iwork1 [nI]
    }

    // memory space for Inext, and Iwork1 of size nI or 2*nI
    GrB_Info info = GB_Work_walloc (iworksize, sizeof (int64_t)) ;
    if (info != GrB_SUCCESS)
    { 
        // out of memory for Work
        GB_wfree ( ) ;
        return (info) ;
    }

    Inext = (int64_t *) GB_thread_local.Work ;

    if (need_Iwork1)
    { 
        // Iwork1 workspace is only needed if the indices I are jumbled,
        // and only for GB_subref_numeric
        Iwork1 = Inext + nI ;            // size nI
    }

    //--------------------------------------------------------------------------
    // ensure Mark is large enough for Mark, of size avlen
    //--------------------------------------------------------------------------

    info = GB_Mark_walloc (avlen) ;
    if (info != GrB_SUCCESS)
    { 
        // out of memory for Mark
        GB_wfree ( ) ;
        return (info) ;
    }

    // ensure flag + nI does not overflow
    Mark = GB_thread_local.Mark ;
    flag = GB_Mark_reset (1, nI) ;

    //--------------------------------------------------------------------------
    // scatter the I indices into buckets
    //--------------------------------------------------------------------------

    // at this point, Mark is clear, so Mark [i] < flag for all i in
    // the range 0 to avlen-1.

    // O(nI) time but this is OK since nI = length of the explicit list I
    for (int64_t inew = nI-1 ; inew >= 0 ; inew--)
    {
        int64_t i = I [inew] ;
        int64_t ihead = (Mark [i] - flag) ;
        if (ihead < 0)
        { 
            // first time i has been seen in the list I
            ihead = -1 ;
        }
        else
        { 
            // i has already been seen in the list I
            nduplicates++ ;
        }
        Mark [i] = inew + flag ;       // (Mark [i] - flag) = inew
        Inext [inew] = ihead ;
    }

    // indices in I are now in buckets.  An index i might appear
    // more than once in the list I.  inew = (Mark [i] - flag) is the
    // first position of i in I (i will be I [inew]), (Mark [i] -
    // flag) is the head of a link list of all places where i appears
    // in I.  inew = Inext [inew] traverses this list, until inew is -1.

    // to traverse all entries in bucket i, do:
    // GB_for_each_entry_in_bucket (inew,i)) { ... }

    #define GB_for_each_entry_in_bucket(inew,i) \
        for (int64_t inew = Mark[i]-flag ; inew >= 0 ; inew = Inext [inew])

    // If Mark [i] < flag, then the ith bucket is empty and i is not in I.
    // Otherise, the first index in bucket i is (Mark [i] - flag).

    #ifndef NDEBUG
    // no part of this code takes O(avlen) time, except this debug test
    for (int64_t i = 0 ; i < avlen ; i++)
    {
        GB_for_each_entry_in_bucket (inew, i)
        {
            ASSERT (inew >= 0 && inew < nI) ;
            ASSERT (i == I [inew]) ;
        }
    }
    #endif

    //--------------------------------------------------------------------------
    // return result
    //--------------------------------------------------------------------------

    *p_Mark = Mark ;
    *p_Inext = Inext ;
    *p_Iwork1 = Iwork1 ;
    *p_nduplicates = nduplicates ;
    *p_flag = flag ;

    return (GrB_SUCCESS) ;
}

