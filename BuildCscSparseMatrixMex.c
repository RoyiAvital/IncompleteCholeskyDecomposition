#include <math.h>
#include <string.h>

#include "mex.h"
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	unsigned int numRows, numCols, numElementsV, numElementsR, numElementsP;
	mwIndex *Ir, *Jc; // Row / Column Vectors of Sparse Matrix
	mwIndex ii;
	double *vV, *vR, *vP, *vData;

	// This function validated the indices and indptr of Numpy / Scipy Sparse match MATLAB's Ir and Jc.
	// See 
	
	// Validate Input
	if( nrhs != 5 || !mxIsNumeric(prhs[0]) || !mxIsNumeric(prhs[1]) || !mxIsNumeric(prhs[2]) ||
		(mxGetNumberOfElements(prhs[0]) != 1) || (mxGetNumberOfElements(prhs[1]) != 1) ||
		(mxGetN(prhs[2]) != 1) || (mxGetN(prhs[3]) != 1) || (mxGetN(prhs[4]) != 1))
	{
		mexErrMsgIdAndTxt("BuildCscSparseMatrix:BuildCscSparseMatrix:nrhs", "Input must be 2 numeric real double scalars and 3 numeric real double column vectors");
	}
	
	// Verify Output
	if( nlhs != 1 ) 
	{
		mexErrMsgIdAndTxt("BuildCscSparseMatrix:BuildCscSparseMatrix:nlhs", "Output must be single sparse matrix");
	}

	mexPrintf("Validated Data\n");
	
	// As in MATLAB values are Doubles
	numRows = (unsigned int)mxGetScalar(prhs[0]);
	numCols = (unsigned int)mxGetScalar(prhs[1]);

	numElementsV = (unsigned int)mxGetM(prhs[2]);
	vV = (double*)mxGetData(prhs[2]);
	numElementsR = (unsigned int)mxGetM(prhs[3]);
	vR = (double*)mxGetData(prhs[3]);
	numElementsP = (unsigned int)mxGetM(prhs[4]);
	vP = (double*)mxGetData(prhs[4]);

	// Verify Dimensions
	if ( (numElementsV != numElementsR) || (numElementsP != numCols + 1) )
	{
		mexErrMsgIdAndTxt("BuildCscSparseMatrix:BuildCscSparseMatrix:nrhs", "The length of vV and vR must match, the length of vP must be numCols + 1");
	}
	
	//// Allocating the output Sparse Matrix
	plhs[0] = mxCreateSparse(numElementsV, numRows, numCols, mxREAL);
	
	vData = (double *)mxGetData(plhs[0]);
	
	 // Allocated on creation of the Sparse Matrix, Length as the number of non zero elements
	Ir = mxGetIr(plhs[0]);
	// Allocated on creation of the Sparse Matrix, Length as the number of columns + 1
	Jc = mxGetJc(plhs[0]);
	
	for (ii = 0; ii < numElementsV; ii++)
	{
		vData[ii]	= vV[ii];
		Ir[ii]		= (mwIndex)vR[ii];
	}

	for (ii = 0; ii < numElementsP; ii++)
	{
		Jc[ii] = (mwIndex)vP[ii];
	}
	
}

