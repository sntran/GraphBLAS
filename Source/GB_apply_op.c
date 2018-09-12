//------------------------------------------------------------------------------
// GB_apply_op:  apply a unary operator to an array
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2018, All Rights Reserved.
// http://suitesparse.com   See GraphBLAS/Doc/License.txt for license.

//------------------------------------------------------------------------------

// Cx = op ((xtype) Ax)

// Compare with GB_transpose_op.c

#include "GB.h"

void GB_apply_op            // apply a unary operator, Cx = op ((xtype) Ax)
(
    void *Cx,               // output array, of type op->ztype
    const GrB_UnaryOp op,   // operator to apply
    const void *Ax,         // input array, of type atype
    const GrB_Type atype,   // type of Ax
    const int64_t anz       // size of Ax and Cx
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    ASSERT (Cx != NULL) ;
    ASSERT (Ax != NULL) ;
    ASSERT (anz >= 0) ;
    ASSERT_OK (GB_check (atype, "atype for GB_apply_op", D0)) ;
    ASSERT_OK (GB_check (op, "op for GB_apply_op", D0)) ;

    //--------------------------------------------------------------------------
    // define the worker for the switch factory
    //--------------------------------------------------------------------------

    // For built-in types only, thus xtype == ztype, but atype can differ
    #define WORKER(ztype,atype)                                 \
    {                                                           \
        ztype *cx = (ztype *) Cx ;                              \
        atype *ax = (atype *) Ax ;                              \
        for (int64_t p = 0 ; p < anz ; p++)                     \
        {                                                       \
            /* x = (ztype) ax [p], type casting */              \
            ztype x ;                                           \
            CAST (x, ax [p]) ;                                  \
            /* apply the unary operator */                      \
            cx [p] = OP (x) ;                                   \
        }                                                       \
        return ;                                                \
    }

    //--------------------------------------------------------------------------
    // launch the switch factory
    //--------------------------------------------------------------------------

    // If GB_COMPACT is defined, the switch factory is disabled and all
    // work is done by the generic worker.  The compiled code will be more
    // compact, but 3 to 4 times slower.

    #ifndef GBCOMPACT

        // switch factory for two types, controlled by code1 and code2
        GB_Type_code code1 = op->ztype->code ;      // defines ztype
        GB_Type_code code2 = atype->code ;          // defines atype

        ASSERT (code1 <= GB_UDT_code) ;
        ASSERT (code2 <= GB_UDT_code) ;

        // BOP(x) is for boolean x, IOP(x) for integer (int* and uint*),
        // and FOP(x) is for floating-point

        switch (op->opcode)
        {

            case GB_ONE_opcode :       // z = 1

                #define BOP(x) true
                #define IOP(x) 1
                #define FOP(x) 1
                #include "GB_2type_template.c"

            case GB_IDENTITY_opcode :  // z = x

                // Do not create workers when the two codes are the same,
                // C is a pure shallow copy of A, and the function has already
                // returned the result C.
                #define NSAME
                #define BOP(x) x
                #define IOP(x) x
                #define FOP(x) x
                #include "GB_2type_template.c"

            case GB_AINV_opcode :      // z = -x

                #define BOP(x)  x
                #define IOP(x) -x
                #define FOP(x) -x
                #include "GB_2type_template.c"

            case GB_ABS_opcode :       // z = abs(x)

                #define BOP(x)  x
                #define IOP(x)  IABS(x)
                #define FOP(x)  FABS(x)
                #include "GB_2type_template.c"

            case GB_MINV_opcode :      // z = 1/x

                // see Source/GB.h discussion on boolean and integer division
                #define BOP(x) true
                #define IOP(x) IMINV(x)
                #define FOP(x) 1./x
                #include "GB_2type_template.c"

            case GB_LNOT_opcode :      // z = ! (x != 0)

                #define BOP(x) !x
                #define IOP(x) (!(x != 0))
                #define FOP(x) (!(x != 0))
                #include "GB_2type_template.c"

            default: ;
        }

    #endif

    #undef WORKER

    // If the switch factory has no worker for the opcode or type, then it
    // falls through to the generic worker below.

    //--------------------------------------------------------------------------
    // generic worker:  apply an operator, with optional typecasting
    //--------------------------------------------------------------------------

    // The generic worker can handle any operator and any type, and it does all
    // required typecasting.  Thus the switch factory can be disabled, and the
    // code will more compact and still work.  It will just be slower.

    int64_t asize = atype->size ;
    int64_t zsize = op->ztype->size ;
    GB_cast_function
        cast_A_to_X = GB_cast_factory (op->xtype->code, atype->code) ;
    GB_unary_function fop = op->function ;

    // scalar workspace
    char xwork [op->xtype->size] ;

    for (int64_t p = 0 ; p < anz ; p++)
    { 
        // xwork = (xtype) Ax [p]
        cast_A_to_X (xwork, Ax +(p*asize), asize) ;
        // Cx [p] = fop (xwork)
        fop (Cx +(p*zsize), xwork) ;
    }
}

