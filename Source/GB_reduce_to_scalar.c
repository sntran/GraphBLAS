//------------------------------------------------------------------------------
// GB_reduce_to_scalar: reduce a matrix to a scalar
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017, All Rights Reserved.
// http://suitesparse.com   See GraphBLAS/Doc/License.txt for license.

//------------------------------------------------------------------------------

// c = accum (c, reduce_to_scalar(A)), reduce entries in a matrix
// to a scalar.  Not user-callable.  Does the work for GrB_*_reduce_TYPE,
// both matrix and vector.

#include "GB.h"

GrB_Info GB_reduce_to_scalar    // twork = reduce_to_scalar (A)
(
    void *c,                    // result scalar
    const GrB_Type ctype,       // the type of scalar, c
    const GrB_BinaryOp accum,   // for c = accum(c,twork)
    const GrB_Monoid reduce,    // monoid to do the reduction
    const GrB_Matrix A          // matrix to reduce
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    // delete any lingering zombies and assemble any pending tuples
    // (required by Table 2.4 of the spec)
    APPLY_PENDING_UPDATES (A) ;
    RETURN_IF_NULL_OR_UNINITIALIZED (reduce) ;
    RETURN_IF_UNINITIALIZED (accum) ;
    RETURN_IF_NULL (c) ;

    ASSERT_OK (GB_check (ctype, "type of scalar c", 0)) ;
    ASSERT_OK (GB_check (reduce, "reduce for reduce_to_scalar", 0)) ;
    ASSERT_OK_OR_NULL (GB_check (accum, "accum for reduce_to_scalar", 0)) ;
    ASSERT_OK (GB_check (A, "A for reduce_to_scalar", 0)) ;

    // check domains and dimensions for c = accum (c,twork)
    GrB_Type ztype = reduce->op->ztype ;
    GrB_Info info = GB_compatible (ctype, NULL, NULL, accum, ztype) ;
    if (info != GrB_SUCCESS)
    {
        return (info) ;
    }

    // twork = reduce (twork,A) must be compatible
    if (!GB_Type_compatible (A->type, ztype))
    {
        return (ERROR (GrB_DOMAIN_MISMATCH, (LOG,
            "incompatible type for reduction operator z=%s(x,y):\n"
            "input of type [%s]\n"
            "cannot be typecast to reduction operator of type [%s]",
            reduce->op->name, A->type->name, reduce->op->ztype->name))) ;
    }

    //--------------------------------------------------------------------------
    // scalar workspace
    //--------------------------------------------------------------------------

    int64_t asize = A->type->size ;
    int64_t anz = NNZ (A) ;
    const void *Ax = A->x ;

    int64_t zsize = ztype->size ;

    char awork [zsize] ;
    char wwork [zsize] ;
    char twork [zsize] ;

    //--------------------------------------------------------------------------
    // twork = reduce_to_scalar (A)
    //--------------------------------------------------------------------------

    // twork = 0
    memcpy (twork, reduce->identity, zsize) ;

    // reduce all the entries in the matrix

    if (A->type == ztype)
    {

        //----------------------------------------------------------------------
        // sum up the entries; no casting needed
        //----------------------------------------------------------------------

        // There are 44 common cases of this function for built-in types and
        // operators.  Four associative operators: min, max, plus, and times
        // with 10 types (int*, uint*, float, and double), and three logical
        // operators (or, and, xor, eq) with a boolean type of C.  All 44 are
        // hard-coded below via a switch factory.  If the case is not handled
        // by the switch factory, 'done' remains false.

        bool done = false ;

        // define the worker for the switch factory
        #define WORKER(type)                                                \
        {                                                                   \
            const type *ax = (type *) Ax ;                                  \
            type s ;                                                        \
            memcpy (&s, twork, zsize) ;                                     \
            for (int64_t p = 0 ; p < anz ; p++)                             \
            {                                                               \
                /* s "+=" ax [p] */                                         \
                ADD (s, ax [p]) ;                                           \
            }                                                               \
            memcpy (twork, &s, zsize) ;                                     \
            done = true ;                                                   \
        }

        //----------------------------------------------------------------------
        // launch the switch factory
        //----------------------------------------------------------------------

        // If GBCOMPACT is defined, the switch factory is disabled and all
        // work is done by the generic worker.  The compiled code will be more
        // compact, but 3 to 4 times slower.

        #ifndef GBCOMPACT

            // controlled by opcode and typecode
            GB_Opcode opcode = reduce->op->opcode ;
            GB_Type_code typecode = A->type->code ;
            ASSERT (typecode <= GB_UDT_code) ;
            #include "GB_assoc_template.c"

        #endif

        #undef WORKER

        //----------------------------------------------------------------------
        // generic worker: sum up the entries, no typecasting
        //----------------------------------------------------------------------

        if (!done)
        {

            GB_binary_function freduce = reduce->op->function ;

            // the switch factory didn't handle this case
            for (int64_t p = 0 ; p < anz ; p++)
            {
                // wwork = twork
                memcpy (wwork, twork, zsize) ;
                // twork = wwork "+" Ax [p]
                freduce (twork, wwork, Ax +(p*asize)) ;
            }
        }

    }
    else
    {

        //----------------------------------------------------------------------
        // generic worker: sum up the entries, with typecasting
        //----------------------------------------------------------------------

        GB_binary_function freduce = reduce->op->function ;
        GB_cast_function
            cast_A_to_Z = GB_cast_factory (ztype->code, A->type->code) ;

        for (int64_t p = 0 ; p < anz ; p++)
        {
            // awork = (ztype) Ax [p]
            cast_A_to_Z (awork, Ax +(p*asize), zsize) ;
            // wwork = twork
            memcpy (wwork, twork, zsize) ;
            // twork = wwork "+" awork
            freduce (twork, wwork, awork) ;
        }
    }

    //--------------------------------------------------------------------------
    // c = twork or c = accum (c,twork)
    //--------------------------------------------------------------------------

    // This operation does not use GB_accum_mask, since c and twork are
    // scalars, not matrices.  There is no scalar mask.

    if (accum == NULL)
    {
        // c = (ctype) twork
        GB_cast_function
            cast_Z_to_C = GB_cast_factory (ctype->code, ztype->code) ;
        cast_Z_to_C (c, twork, ctype->size) ;
    }
    else
    {
        GB_binary_function faccum = accum->function ;

        GB_cast_function cast_C_to_xaccum, cast_Z_to_yaccum, cast_zaccum_to_C ;
        cast_C_to_xaccum = GB_cast_factory (accum->xtype->code, ctype->code) ;
        cast_Z_to_yaccum = GB_cast_factory (accum->ytype->code, ztype->code) ;
        cast_zaccum_to_C = GB_cast_factory (ctype->code, accum->ztype->code) ;

        // scalar workspace
        char xaccum [accum->xtype->size] ;
        char yaccum [accum->ytype->size] ;
        char zaccum [accum->ztype->size] ;

        // xaccum = (accum->xtype) c
        cast_C_to_xaccum (xaccum, c, ctype->size) ;

        // yaccum = (accum->ytype) twork
        cast_Z_to_yaccum (yaccum, twork, zsize) ;

        // zaccum = xaccum "+" yaccum
        faccum (zaccum, xaccum, yaccum) ;

        // c = (ctype) zaccum
        cast_zaccum_to_C (c, zaccum, ctype->size) ;
    }

    return (REPORT_SUCCESS) ;
}

