//------------------------------------------------------------------------------
// GrB_init: initialize GraphBLAS
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2018, All Rights Reserved.
// http://suitesparse.com   See GraphBLAS/Doc/License.txt for license.

//------------------------------------------------------------------------------

// GrB_init must called before any other GraphBLAS operation.  GrB_finalize
// must be called as the last GraphBLAS operation.

// GrB_init defines the mode that GraphBLAS will use:  blocking or
// non-blocking.  With blocking mode, all operations finish before returning to
// the user application.  With non-blocking mode, operations can be left
// pending, and are computed only when needed.

// The GrB_wait function forces all pending operations to complete.  Blocking
// mode is as if the GrB_wait operation is called whenever a GraphBLAS
// operation returns to the user.

// The non-blocking mode can have side effects if user-defined functions have
// side effects or if they rely on global variables, which are not under the
// control of GraphBLAS.  Suppose the user creates a user-defined operator that
// accesses a global variable.  That operator is then used in a GraphBLAS
// operation, which is left pending.  If the user then changes the global
// variable, the pending operations will be eventually computed with this
// different value.

// Worse yet, a user-defined operator can be freed before it is needed to
// finish a pending operation.  To avoid this, call GrB_wait before modifying
// any global variables relied upon by user-defined operators and before
// freeing any user-defined types, operators, monoids, or semirings.

// This version of SuiteSparse:GraphBLAS does not actually require a call to
// GrB_init.  All required global variables are statically initialized to their
// proper values.  However, for best practice, GrB_init should be called prior
// to any other call to GraphBLAS functions.

#include "GB.h"

//------------------------------------------------------------------------------
// All thread local storage is held in a single struct, initialized here
//------------------------------------------------------------------------------

_Thread_local GB_thread_local_struct GB_thread_local =
{
    // error status
    .info = GrB_SUCCESS,        // last info status of user-callable routines
    .row = 0,                   // last row index searched for
    .col = 0,                   // last column index searched for
    .is_matrix = 0,             // last search matrix (true) or vector (false)
    .where = "",                // last user-callable function called
    .file = "",                 // file where error occurred
    .line = 0,                  // line number where error occurred
    .details = "",              // details of the error
    .report = "",               // report created by GrB_error

    // thread private workspace
    .Mark = NULL,               // initialized space
    .Mark_flag = 1,             // current watermark in Mark [...]
    .Mark_size = 0,             // current size of Mark array
    .Work = NULL,               // uninitialized space
    .Work_size = 0,             // current size of Work array
    .Flag = NULL,               // initialized space
    .Flag_size = 0,             // size of Flag array

    // random seed for each thread
    .seed = 1,

    // last method used in C=A*B
    .AxB_method = GxB_DEFAULT
} ;

//------------------------------------------------------------------------------
// All Global storage is declared and initialized here
//------------------------------------------------------------------------------

// If the user creates threads that work on GraphBLAS matrices, then all of
// those threads must share the same matrix queue, and the same mode.

GB_Global_struct GB_Global =
{

    // queued matrices with work to do
    .queue_head = NULL,         // pointer to first queued matrix

    // GraphBLAS mode
    .mode = GrB_NONBLOCKING,    // default is nonblocking

    // initialization flag
    .GrB_init_called = false,   // GrB_init has not yet been called

    // default format
    .hyper_ratio = GB_HYPER_DEFAULT,
    .is_csc = (GB_FORMAT_DEFAULT != GxB_BY_ROW),    // default is GxB_BY_ROW

    // malloc tracking, for testing, statistics, and debugging only
    .nmalloc = 0,               // memory block counter
    .malloc_debug = false,      // do not test memory handling
    .malloc_debug_count = 0,    // counter for testing memory handling
    .inuse = 0,                 // memory space current in use
    .maxused = 0                // high water memory usage
} ;

//------------------------------------------------------------------------------
// GrB_init
//------------------------------------------------------------------------------

// If GraphBLAS is used by multiple user threads, only one can call GrB_init.

GrB_Info GrB_init           // start up GraphBLAS
(
    const GrB_Mode mode     // blocking or non-blocking mode
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    GB_WHERE ("GrB_init (mode)") ;

    if (! (mode == GrB_BLOCKING || mode == GrB_NONBLOCKING))
    { 
        return (GB_ERROR (GrB_INVALID_VALUE, (GB_LOG,
            "Unknown mode: %d; must be %d (GrB_NONBLOCKING) or %d"
            " (GrB_BLOCKING)", mode, GrB_NONBLOCKING, GrB_BLOCKING))) ;
    }

    //--------------------------------------------------------------------------
    // initialize per-thread workspace
    //--------------------------------------------------------------------------

    // GrB_init should be called at least once by the user application.  It is
    // safe for any thread to skip this, since the settings are properly
    // defined above.

    // per-thread error status
    GB_thread_local.info = GrB_SUCCESS ;
    GB_thread_local.row = 0 ;
    GB_thread_local.col = 0 ;
    GB_thread_local.is_matrix = 0 ;
    GB_thread_local.file = __FILE__ ;
    GB_thread_local.line = __LINE__ ;
    GB_thread_local.details [0] = '\0' ;
    GB_thread_local.report  [0] = '\0' ;

    // workspace
    GB_thread_local.Mark = NULL ;       // initialized space
    GB_thread_local.Mark_flag = 1 ;     // current watermark in Mark [...]
    GB_thread_local.Mark_size = 0 ;     // current size of Mark array
    GB_thread_local.Work = NULL ;       // uninitialized space
    GB_thread_local.Work_size = 0 ;     // current size of Work array
    GB_thread_local.Flag = NULL ;       // initialized space
    GB_thread_local.Flag_size = 0 ;     // current size of Flag array

    // random seed
    GB_thread_local.seed = 1 ;

    // last method used in C=A*B
    GB_thread_local.AxB_method = GxB_DEFAULT ;

    //--------------------------------------------------------------------------
    // initialize GraphBLAS global status
    //--------------------------------------------------------------------------

    // Only one thread should initialize these settings.  If multiple threads
    // call GrB_init, only the first thread does this work.

    bool I_was_first = false ;

    // queue of matrices for nonblocking mode and set the mode.  The queue
    // must be protected and can be initialized only once by any thread.
    #pragma omp critical (GB_queue)
    {
        if (!GB_Global.GrB_init_called)
        { 
            I_was_first = true ;

            // clear the queue
            GB_Global.queue_head = NULL ;

            // set the mode: blocking or nonblocking
            GB_Global.mode = mode ;             // default is non-blocking

            // first thread has called GrB_init
            GB_Global.GrB_init_called = true ;
        }
    }

    if (! I_was_first)
    { 
        return (GB_ERROR (GrB_INVALID_VALUE, (GB_LOG,
            "GrB_init must not be called twice"))) ;
    }

    // set the default hypersparsity ratio and CSR/CSC format;  any thread
    // can do this later as well, so there is no race condition danger.
    GB_global_option_set (true, GB_HYPER_DEFAULT, true, GB_FORMAT_DEFAULT) ;

    // malloc tracking.  This is only for statistics and development.
    // FUTURE: these statistics could be per-thread.
    #pragma omp critical (GB_memory)
    {
        GB_Global.nmalloc = 0 ;
        GB_Global.malloc_debug = false ;
        GB_Global.malloc_debug_count = 0 ;
        GB_Global.inuse = 0 ;
        GB_Global.maxused = 0 ;
    }
    return (GB_REPORT_SUCCESS) ;
}

