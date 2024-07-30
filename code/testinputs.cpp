/*
 * testinputs.cpp
 */

#include "mex.h"


void mexFunction(int  nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	// This function will error if number of inputs its not 3 or 4
	// This function will error if number of outputs is more than 1

	    // Check inputs:
	    if (nrhs < 3 || nrhs > 4) {
	        mexErrMsgIdAndTxt("Testinputs:ErrorIdIn",
	            "Invalid number of inputs to MEX file.");
	    }

	    // Check outputs:
	    if (nlhs > 1) {
	        mexErrMsgIdAndTxt("Testinputs:ErrorIdOut",
	                "Invalid number of outputs to MEX file.");
	    }

}

