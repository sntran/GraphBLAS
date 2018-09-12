//------------------------------------------------------------------------------
// GB_matvec_check: print a GraphBLAS matrix and check if it is valid
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2018, All Rights Reserved.
// http://suitesparse.com   See GraphBLAS/Doc/License.txt for license.

//------------------------------------------------------------------------------

#include "GB.h"

GrB_Info GB_matvec_check    // check a GraphBLAS matrix or vector
(
    const GrB_Matrix A,     // GraphBLAS matrix to print and check
    const char *name,       // name of the matrix, optional
    int pr,                 // 0: print nothing, 1: print header and errors,
                            // 2: print brief, 3: print all
                            // if negative, ignore queue conditions
                            // and use FLIP(pr) for diagnostic printing.
    const char *kind        // "matrix" or "vector"
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    bool ignore_queue = false ;
    if (pr < 0)
    { 
        // -2: print nothing (pr = 0)
        // -3: print header  (pr = 1)
        // -4: print brief   (pr = 2)
        // -5: print all     (pr = 3)
        pr = FLIP (pr) ;
        ignore_queue = true ;
    }

    if (pr > 0) printf ("\nGraphBLAS %s: %s ", kind, NAME) ;

    if (A == NULL)
    { 
        // GrB_error status not modified since this may be an optional argument
        if (pr > 0) printf ("NULL\n") ;
        return (GrB_NULL_POINTER) ;
    }

    //--------------------------------------------------------------------------
    // check the object
    //--------------------------------------------------------------------------

    CHECK_MAGIC (A, kind) ;
    ASSERT (A->magic == MAGIC) ;    // A is now a valid initialized object

    if (pr > 0)
    { 
        printf ("\nnrows: "GBd" ncols: "GBd" max # entries: "GBd"\n",
            NROWS (A), NCOLS (A), A->nzmax) ;
        printf ("format: %s %s",
            A->is_hyper ? "hypersparse" : "standard",
            A->is_csc ?   "CSC" : "CSR") ;
        printf (" vlen: "GBd" nvec_nonempty: "GBd" nvec: "GBd" plen: "
            GBd " vdim: "GBd"\n",
            A->vlen, A->nvec_nonempty, A->nvec, A->plen, A->vdim) ;
        printf ("hyper_ratio %g\n", A->hyper_ratio) ;
    }

    if (A->vlen < 0 || A->vlen > GB_INDEX_MAX ||
        A->vdim < 0 || A->vdim > GB_INDEX_MAX ||
        A->nzmax < 0 || A->nzmax > GB_INDEX_MAX)
    { 
        if (pr > 0) printf ("invalid %s dimensions\n", kind) ;
        return (ERROR (GrB_INVALID_OBJECT, (LOG,
            "%s invalid : nrows, ncols, or nzmax out of range: [%s]",
            kind, NAME))) ;
    }

    if (A->is_hyper)
    {
        if (! (A->nvec >= 0 && A->nvec <= A->plen && A->plen <= A->vdim &&
               A->nvec == A->nvec_nonempty))
        { 
            if (pr > 0) printf ("invalid %s hypersparse structure\n", kind) ;
            return (ERROR (GrB_INVALID_OBJECT, (LOG,
                "%s invalid hypersparse structure [%s]", kind, NAME))) ;
        }
    }
    else
    {
        if (! (A->nvec == A->plen && A->plen == A->vdim))
        { 
            if (pr > 0) printf ("invalid %s standard structure\n", kind) ;
            return (ERROR (GrB_INVALID_OBJECT, (LOG,
                "%s invalid structure [%s]", kind, NAME))) ;
        }
    }

    // a matrix contains 1 to 8 different malloc'd blocks
    int64_t nallocs = 1 +                       // header
        (A->h != NULL && !A->h_shallow) +       // A->h, if not shallow
        (A->p != NULL && !A->p_shallow) +       // A->p, if not shallow
        (A->i != NULL && !A->i_shallow) +       // A->i, if not shallow
        (A->x != NULL && !A->x_shallow) +       // A->x, if not shallow
        (A->i_pending != NULL) +                // A->i_pending if tuples
        (A->j_pending != NULL) +                // A->j_pending if tuples
        (A->s_pending != NULL) ;                // A->s_pending if tuples

    #ifdef DEVELOPER
    if (pr > 1) printf ("number of memory blocks: "GBd"\n", nallocs) ;
    #endif

    GrB_Info info = GB_check (A->type, "", pr) ;
    if (info != GrB_SUCCESS)
    { 
        if (pr > 0) printf ("%s has an invalid type\n", kind) ;
        return (ERROR (GrB_INVALID_OBJECT, (LOG,
            "%s has an invalid type: [%s]", kind, NAME))) ;
    }

    #ifdef DEVELOPER
    if (pr > 1) printf ("->h: %p shallow: %d\n", A->h, A->h_shallow) ;
    if (pr > 1) printf ("->p: %p shallow: %d\n", A->p, A->p_shallow) ;
    if (pr > 1) printf ("->i: %p shallow: %d\n", A->i, A->i_shallow) ;
    if (pr > 1) printf ("->x: %p shallow: %d\n", A->x, A->x_shallow) ;
    #endif

    if (A->p == NULL)
    { 
        if (pr > 0) printf ("->p is NULL, invalid %s\n", kind) ;
        return (ERROR (GrB_INVALID_OBJECT, (LOG,
            "%s contains a NULL A->p pointer: [%s]", kind, NAME))) ;
    }

    if (A->is_hyper)
    {
        if (A->h == NULL)
        { 
            if (pr > 0) printf ("->h is NULL, invalid hypersparse %s\n",
                kind) ;
            return (ERROR (GrB_INVALID_OBJECT, (LOG,
                "hypersparse %s contains a NULL A->h pointer: [%s]",
                kind, NAME))) ;
        }
    }
    else
    {
        if (A->h != NULL)
        { 
            if (pr > 0) printf ("->h is not NULL, invalid non-hypersparse %s\n",
                kind) ;
            return (ERROR (GrB_INVALID_OBJECT, (LOG,
                "non-hypersparse %s contains a non-NULL A->h pointer: [%s]",
                kind, NAME))) ;
        }
    }

    bool A_empty = (A->nzmax == 0) ;

    if (A_empty)
    {
        // A->x and A->i pointers must be NULL and shallow must be false
        if (A->i != NULL || A->i_shallow || A->x_shallow)
        { 
            if (pr > 0) printf ("invalid empty %s\n", kind) ;
            return (ERROR (GrB_INVALID_OBJECT, (LOG,
                "%s is an invalid empty object: [%s]", kind, NAME))) ;
        }

        // check the vector pointers
        for (int64_t j = 0 ; j <= A->nvec ; j++)
        {
            if (A->p [j] != 0)
            { 
                if (pr > 0) printf ("->p ["GBd"] = "GBd" invalid\n", j,A->p[j]);
                return (ERROR (GrB_INVALID_OBJECT, (LOG,
                    "%s ->p ["GBd"] = "GBd" invalid: [%s]",
                    kind, j, A->p[j], NAME))) ;
            }
        }
        if (pr > 0) printf ("empty\n") ;
    }

    if (!A_empty && A->i == NULL)
    { 
        if (pr > 0) printf ("->i is NULL, invalid %s\n", kind) ;
        return (ERROR (GrB_INVALID_OBJECT, (LOG,
            "%s contains a NULL A->i pointer: [%s]", kind, NAME))) ;
    }

    //--------------------------------------------------------------------------
    // check the vector pointers
    //--------------------------------------------------------------------------

    if (A->p [0] != 0)
    { 
        if (pr > 0) printf ("->p [0] = "GBd" invalid\n", A->p [0]) ;
        return (ERROR (GrB_INVALID_OBJECT, (LOG,
            "%s A->p [0] = "GBd" invalid: [%s]", kind, A->p [0], NAME))) ;
    }

    for (int64_t j = 0 ; j < A->nvec ; j++)
    {
        if (A->p [j+1] < A->p [j] || A->p [j+1] > A->nzmax)
        { 
            if (pr > 0) printf ("->p ["GBd"] = "GBd" invalid\n",
                j+1, A->p [j+1]) ;
            return (ERROR (GrB_INVALID_OBJECT, (LOG,
                "%s A->p ["GBd"] = "GBd" invalid: [%s]",
                kind, j+1, A->p [j+1], NAME))) ;
        }
    }

    if (A->is_hyper)
    {
        int64_t jlast = -1 ;
        for (int64_t k = 0 ; k < A->nvec ; k++)
        {
            int64_t j = A->h [k] ;
            if (jlast >= j || j < 0 || j >= A->vdim)
            { 
                if (pr > 0) printf ("->h ["GBd"] = "GBd" invalid\n",
                    k, A->h [k]) ;
                return (ERROR (GrB_INVALID_OBJECT, (LOG,
                    "%s A->h ["GBd"] = "GBd" invalid: [%s]",
                    kind, k, A->h [k], NAME))) ;
            }
            jlast = j ;
        }
    }

    int64_t anz = NNZ (A) ;
    if (pr > 0) printf ("number of entries: "GBd" ", anz) ;

    if (pr > 0) printf ("\n") ;

    //--------------------------------------------------------------------------
    // report the number of pending tuples and number of zombies
    //--------------------------------------------------------------------------

    if (A->n_pending != 0 || A->nzombies != 0)
    { 
        if (pr > 0) printf ("pending tuples: "GBd" max pending: "GBd
            " zombies: "GBd"\n", A->n_pending, A->max_n_pending, A->nzombies) ;
    }

    if (A->nzombies < 0 || A->nzombies > anz)
    { 
        if (pr > 0) printf ("invalid number of zombies: "GBd" "
            "must be >= 0 and <= # entries ("GBd")\n", A->nzombies, anz) ;
        return (ERROR (GrB_INVALID_OBJECT, (LOG,
            "%s invalid number of zombies: "GBd"\n"
            "must be >= 0 and <= # entries ("GBd") [%s]",
            kind, A->nzombies, anz, NAME))) ;
    }

    //--------------------------------------------------------------------------
    // check and print the row indices and numerical values
    //--------------------------------------------------------------------------

    #define NBRIEF 10
    #define NZBRIEF 30

    bool jumbled = false ;
    int64_t nzombies = 0 ;
    int64_t jcount = 0 ;

    for_each_vector (A)
    {
        int64_t ilast = -1 ;
        for_each_entry (j, p, pend)
        {
            bool prcol = ((pr > 1 && jcount < NBRIEF) || pr > 2) ;
            if (ilast == -1)
            {
                // print the header for vector j
                if (prcol)
                { 
                    printf ("%s: "GBd" : "GBd" entries ["GBd":"GBd"]\n",
                        A->is_csc ? "column" : "row", j, pend - p, p, pend-1) ;
                }
                else if (pr == 2 && jcount == NBRIEF)
                { 
                    printf ("...\n") ;
                }
                jcount++ ;      // count # of vectors printed so far
            }
            int64_t i = A->i [p] ;
            bool is_zombie = IS_ZOMBIE (i) ;
            i = UNFLIP (i) ;
            if (is_zombie) nzombies++ ;
            if (prcol)
            { 
                if ((pr > 1 && p < NZBRIEF) || pr > 2)
                { 
                    printf ("    %s "GBd": ", A->is_csc ? "row" : "column", i) ;
                }
                else if (pr == 2 && (ilast == -1 || p == NZBRIEF))
                { 
                    printf ("    ...\n") ;
                }
            }
            int64_t row = A->is_csc ? i : j ;
            int64_t col = A->is_csc ? j : i ;
            if (i < 0 || i >= A->vlen)
            { 
                if (pr > 0) printf ("index ("GBd","GBd") out of range\n",
                    row, col) ;
                return (ERROR (GrB_INVALID_OBJECT, (LOG,
                    "%s index ("GBd","GBd") out of range: [%s]",
                    kind, row, col, NAME))) ;
            }

            // print the value
            bool print_value = (prcol && ((pr > 1 && p < NZBRIEF) || pr > 2)) ;
            if (print_value)
            { 
                if (is_zombie)
                { 
                    printf ("zombie") ;
                }
                else if (A->x != NULL)
                { 
                    GB_Entry_print (A->type, A->x +(p * A->type->size)) ;
                }
            }

            if (i <= ilast)
            { 
                // indices unsorted, or duplicates present
                if (pr > 0) printf (" index ("GBd","GBd") jumbled", row, col) ;
                jumbled = true ;
                print_value = (pr > 0) ;
            }

            if (print_value)
            { 
                printf ("\n") ;
            }
            ilast = i ;
        }
    }

    //--------------------------------------------------------------------------
    // check the zombie count
    //--------------------------------------------------------------------------

    if (nzombies != A->nzombies)
    { 
        if (pr > 0) printf ("invalid zombie count: "GBd" exist but"
            " A->nzombies = "GBd"\n", nzombies, A->nzombies) ;
        return (ERROR (GrB_INVALID_OBJECT, (LOG,
            "%s invalid zombie count: "GBd" exist but A->nzombies = "GBd" "
            "[%s]", kind, nzombies, A->nzombies, NAME))) ;
    }

    //--------------------------------------------------------------------------
    // check and print the pending tuples
    //--------------------------------------------------------------------------

    if (A->n_pending < 0 || A->n_pending > A->max_n_pending ||
        A->max_n_pending < 0)
    { 
        if (pr > 0) printf ("invalid pending count\n") ;
        return (ERROR (GrB_INVALID_OBJECT, (LOG,
            "%s invalid pending tuple count: pending "GBd" max "GBd": [%s]",
            kind, A->n_pending, A->max_n_pending, NAME))) ;
    }

    #ifdef DEVELOPER
    if (pr > 1) printf ("A %p\n", A) ;
    if (pr > 1) printf ("->i_pending %p\n", A->i_pending) ;
    if (pr > 1) printf ("->j_pending %p\n", A->j_pending) ;
    if (pr > 1) printf ("->s_pending %p\n", A->s_pending) ;
    #endif

    if (A->n_pending == 0)
    {

        //---------------------------------------------------------------------
        // A has no pending tuples
        //---------------------------------------------------------------------

        // no tuples; arrays must be NULL
        if (A->i_pending != NULL || A->s_pending != NULL ||
            A->j_pending != NULL || A->max_n_pending != 0)
        { 
            if (pr > 0) printf ("invalid pending tuples\n") ;
            return (ERROR (GrB_INVALID_OBJECT, (LOG,
                "%s invalid pending tuples: [%s]", kind, NAME))) ;
        }

    }
    else
    {

        //---------------------------------------------------------------------
        // A has pending tuples
        //---------------------------------------------------------------------

        // matrix has tuples, arrays and type must not be NULL
        if (A->i_pending == NULL || A->s_pending == NULL ||
            (A->vdim > 1 && A->j_pending == NULL))
        { 
            if (pr > 0) printf ("invalid pending tuples\n") ;
            return (ERROR (GrB_INVALID_OBJECT, (LOG,
                "%s invalid pending tuples: [%s]", kind, NAME))) ;
        }

        if (pr > 0) printf ("pending tuples:\n") ;

        info = GB_check (A->type_pending, "", pr) ;
        if (info != GrB_SUCCESS)
        { 
            if (pr > 0) printf ("%s has an invalid type_pending\n", kind) ;
            return (ERROR (GrB_INVALID_OBJECT, (LOG,
                "%s has an invalid type_pending: [%s]", kind, NAME))) ;
        }

        int64_t ilast = -1 ;
        int64_t jlast = -1 ;
        bool sorted = true ;

        for (int64_t k = 0 ; k < A->n_pending ; k++)
        {
            int64_t i = A->i_pending [k] ;
            int64_t j = (A->vdim <= 1) ? 0 : (A->j_pending [k]) ;
            int64_t row = A->is_csc ? i : j ;
            int64_t col = A->is_csc ? j : i ;

            // print the tuple
            if ((pr > 1 && k < NZBRIEF) || pr > 2)
            { 
                printf ("row: "GBd" col: "GBd" ", row, col) ;
                GB_Entry_print (A->type_pending,
                    A->s_pending +(k * A->type_pending->size)) ;
                printf ("\n") ;
            }

            if (i < 0 || i >= A->vlen || j < 0 || j >= A->vdim)
            { 
                if (pr > 0) printf ("tuple ("GBd","GBd") out of range\n",
                    row, col) ;
                return (ERROR (GrB_INVALID_OBJECT, (LOG,
                    "%s tuple index ("GBd","GBd") out of range: [%s]",
                    kind, row, col, NAME))) ;
            }

            sorted = sorted && ((jlast < j) || (jlast == j && ilast <= i)) ;
            ilast = i ;
            jlast = j ;
        }

        if (sorted != A->sorted_pending)
        { 
            printf ("sorted %d sorted_pending %d\n", sorted, A->sorted_pending);
            if (pr > 0) printf ("invalid pending tuples: invalid sort\n") ;
            return (ERROR (GrB_INVALID_OBJECT, (LOG,
                "%s invalid pending tuples: [%s]", kind, NAME))) ;
        }

        if (A->operator_pending == NULL)
        { 
            if (pr > 0) printf ("pending operator: implicit 2nd\n") ;
        }
        else
        {
            info = GB_check (A->operator_pending, "pending operator:", pr) ;
            if (info != GrB_SUCCESS)
            { 
                if (pr > 0) printf ("invalid pending operator\n") ;
                return (ERROR (GrB_INVALID_OBJECT, (LOG,
                    "%s invalid operator: [%s]", kind, NAME))) ;
            }
        }
    }

    //--------------------------------------------------------------------------
    // check the queue
    //--------------------------------------------------------------------------

    if (!ignore_queue)
    {
        GrB_Matrix head, prev, next ;
        bool enqd ;

        GB_queue_check (A, &head, &prev, &next, &enqd) ;

        #ifdef DEVELOPER
        if (pr > 1) printf ("queue head  %p\n", head) ;
        if (pr > 1) printf ("queue prev  %p\n", prev) ;
        if (pr > 1) printf ("queue next  %p\n", next) ;
        if (pr > 1) printf ("is in queue %d\n", enqd) ;
        #endif

        #define IS_NOT_IN_QUEUE(A) (prev == NULL && head != A)
        #define IS_IN_QUEUE(A) (! IS_NOT_IN_QUEUE(A))
        if (enqd != IS_IN_QUEUE (A))
        { 
            if (pr > 0) printf ("queued state inconsistent: [%d] != [%d]\n",
                enqd, IS_IN_QUEUE (A)) ;
            return (ERROR (GrB_INVALID_OBJECT, (LOG,
                "%s queued state inconsistent: [%s], [%d] != [%d]", kind, NAME,
                enqd, IS_IN_QUEUE (A)))) ;
        }
        #undef IS_NOT_IN_QUEUE
        #undef IS_IN_QUEUE

        if (PENDING (A) || ZOMBIES (A))
        {
            if (!enqd)
            { 
                if (pr > 0) printf ("must be in queue but is not there\n") ;
                return (ERROR (GrB_INVALID_OBJECT, (LOG,
                "%s must be in queue but is not there: [%s]", kind, NAME))) ;
            }

            // prev is NULL if and only if A is at the head of the queue
            if ((prev == NULL) != (head == A))
            { 
                if (pr > 0) printf ("invalid queue\n") ;
                return (ERROR (GrB_INVALID_OBJECT, (LOG,
                    "%s invalid queue: [%s]", kind, NAME))) ;
            }
        }
        else
        {
            if (enqd)
            { 
                if (pr > 0) printf ("must not be in queue but is there\n") ;
                return (ERROR (GrB_INVALID_OBJECT, (LOG,
                    "%s must not be in queue but present there: [%s]",
                    kind, NAME))) ;
            }
        }
    }

    if (pr == 3) printf ("\n") ;

    //--------------------------------------------------------------------------
    // return result
    //--------------------------------------------------------------------------

    int64_t actual = GB_nvec_nonempty (A) ;
    if (A->nvec_nonempty != actual)
    { 
        if (pr > 0) printf ("invalid count of non-empty vectors"
            "A->nvec_nonempty = "GBd" actual "GBd"\n",
            A->nvec_nonempty, actual) ;
        return (ERROR (GrB_INVALID_OBJECT, (LOG,
            "%s invalid count of nonempty-vectors [%s]", kind, NAME))) ;
    }

    // Returns GrB_INVALID_OBJECT if a row or column index is out of bounds,
    // since this indicates the object is corrupted.  No valid matrix is ever
    // built with indices out of bounds since the indices are checked when the
    // matrix is built.

    // Returns GrB_INDEX_OUT_OF_BOUNDS if a column has unsorted indices, and
    // perhaps duplicates as well.  For matrices passed back to the user, or
    // obtained from the user, this is an error.  For some matrices internally,
    // the row indices may be jumbled.  These are about to be sorted via qsort
    // or transpose.  In this case, a jumbled matrix is OK.  Duplicates are
    // still an error but this function does not distinguish between the two
    // cases (it would require workspace to do so).  See the
    // ASSERT_OK_OR_JUMBLED macro.

    // do not log error with REPORT_SUCCESS; it may mask an error in the caller
    return (jumbled ? GrB_INDEX_OUT_OF_BOUNDS : GrB_SUCCESS) ;
}

#undef NBRIEF
#undef NZBRIEF

