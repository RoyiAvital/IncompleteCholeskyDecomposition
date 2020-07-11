
#include <stdlib.h>
#include <string.h>

#include "mex.h"
#include "IncompleteCholeskyDecomposition.c"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	mwSize numRows, numCols, numElements, numNzL;
	mwIndex *Ir, *Jc; // Pseudo Row / Column Vectors of Sparse Matrix
	double *vDataL, *vX, *vY;
	unsigned int *vIndicesL, *vIndicesPtrL, ii;
	
	// Validate Input
	if ( nrhs != 2 )
	{
		mexErrMsgIdAndTxt("SolveCholsekyLinearSystem:SolveCholsekyLinearSystem:nrhs", "There must be 2 inputs: A real sparse lower triangle suare matrix of type double and a real column vector of type double");
	}
	// Validating the 1st input
	if( !mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || (mxGetClassID(prhs[0]) != mxDOUBLE_CLASS) || !mxIsSparse(prhs[0]) )
	{
		mexErrMsgIdAndTxt("SolveCholsekyLinearSystem:SolveCholsekyLinearSystem:nrhs", "The first input must be a real sparse positive definite matrix of type double");
	}

	// Validating the 2nd input
	if ( !mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) || (mxGetClassID(prhs[1]) != mxDOUBLE_CLASS) )
	{
		mexErrMsgIdAndTxt("SolveCholsekyLinearSystem:SolveCholsekyLinearSystem:nrhs", "The second input must be a real column vector of type double");
	}

	// Validate Output
	if (nlhs > 1)
	{
		mexErrMsgIdAndTxt("SolveCholsekyLinearSystem:SolveCholsekyLinearSystem:nlhs", "The function works either in place or with a single output");
	}
	
	// Input Objects Dimensions
	numRows		= mxGetM(prhs[0]);
	numCols		= mxGetN(prhs[0]);
	numElements = mxGetM(prhs[1]);

	// Validation Dimensions
	if ( (numRows != numCols) || (mxGetNumberOfElements(prhs[0]) != numRows * numCols) || (mxGetNumberOfElements(prhs[1]) != numCols) || (mxGetN(prhs[1]) != 1) || (numElements != numCols) )
	{
		mexErrMsgIdAndTxt("SolveCholsekyLinearSystem:SolveCholsekyLinearSystem:nrhs", "The dimensions of the 1st input must be NxN and the 2nd input Nx1");
	}
	
	// Input Matrix Data
	vDataL = (double*)mxGetData(prhs[0]);
	// Allocated on creation of the Sparse Matrix, Length as the number of non zero elements
	Ir = mxGetIr(prhs[0]);
	// Allocated on creation of the Sparse Matrix, Length as the number of columns + 1
	Jc = mxGetJc(prhs[0]);

	vY = (double*)mxGetData(prhs[1]);

	if (nlhs == 1)
	{
		plhs[0] = mxCreateDoubleMatrix(numElements, (mwSize)(1), mxREAL);
		vX = (double*)mxGetData(plhs[0]);
		memcpy(vX, vY, numCols * sizeof(double));
	}
	else
	{
		vX = (double*)mxGetData(prhs[1]);
	}

	numNzL			= Jc[numCols]; // The last element of Jx has the actual number of NNZ in the Matrix.
	vIndicesL		= (unsigned int*)malloc(numNzL * sizeof(unsigned int));
	vIndicesPtrL	= (unsigned int*)malloc((numCols + 1) * sizeof(unsigned int));

	for (ii = 0; ii < numNzL; ii++)
	{
		vIndicesL[ii] = (unsigned int)Ir[ii];
	}

	for (ii = 0; ii < numCols + 1; ii++)
	{
		vIndicesPtrL[ii] = (unsigned int)Jc[ii];
	}
	
	_BackSubstitutionL(vX, vDataL, vIndicesL, vIndicesPtrL, (unsigned int)numCols);
	_BackSubstitutionLT(vX, vDataL, vIndicesL, vIndicesPtrL, (unsigned int)numCols);
	
	free(vIndicesL);
	free(vIndicesPtrL);
	
}