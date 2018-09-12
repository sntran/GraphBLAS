//------------------------------------------------------------------------------
// GB_mx_put_time: put the time back to the global MATLAB workspace
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2018, All Rights Reserved.
// http://suitesparse.com   See GraphBLAS/Doc/License.txt for license.

//------------------------------------------------------------------------------

#include "GB_mex.h"

double gbtime = 0, tic [2] = {0,0} ;

void GB_mx_put_time ( )                 // return the time to MATLAB
{

    // create a MATLAB array with the right size
    mxArray * gbresults_matlab = mxCreateNumericMatrix (1, 2,
            mxDOUBLE_CLASS, mxREAL) ;

    // copy the time into the MATLAB array
    double *t = (double *) mxGetData (gbresults_matlab) ;

    t [0] = gbtime ;
    t [1] = GB_thread_local.AxB_method ;

    GB_thread_local.AxB_method = -1 ;
    gbtime = 0 ;

    // put the MATLAB array into the global workspace, overwriting the
    // version that was already there
    mexPutVariable ("global", "GraphBLAS_results", gbresults_matlab) ;
}


