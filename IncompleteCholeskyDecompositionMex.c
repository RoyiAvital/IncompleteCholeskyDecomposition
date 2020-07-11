/*
 * Incomplete Cholesky Decomposition - MEX Wrapper
 * This module implements the Cholesky Decomposition. It uses upper boundary to the number of
 * elements to prevent recurring allocations. It is implemented with the preconditioned conjugate
 * gradient in mind hence the preconditioning step is implemented as well. The decomposed Matrix
 * is given by mA + shiftVal * diag(mA).
 * Input:
 *	- mA                -   Input Positive Definite Sparse Matrix.
 *							The matrix to decompose.
 *                          Structure: Sparse (CSC) Matrix (numRows * numCols).
 *                          Type: 'Double'.
 *                          Range: (-inf, inf).
 *	- discardThr        -   Discard Threshold.
 *							Values in the cholesky decomposition which are smaller than
 *							`discardThr` are zeroed. If `discardThr` = 0 then the full
 *							Cholesky Decomposition is applied (No discarding).
 *                          Structure: Scalar.
 *                          Type: 'Double'.
 *                          Range: [0, inf).
 * - shiftVal         -   Discard Threshold.
 *							Values in the cholesky decomposition which are smaller than
 *							`discardThr` are zeroed.
 *                          Structure: Scalar.
 *                          Type: 'Double'.
 *                          Range: (0, inf).
 * - maxNumNz         -   Discard Threshold.
 *							Values in the cholesky decomposition which are smaller than
 *							`discardThr` are zeroed.
 *                          Structure: Scalar.
 *                          Type: 'Double'.
 *                          Range: (0, inf).
 * Output:
 *	- mL                -   Input Positive Definite Sparse Matrix.
 *							The matrix to decompose.
 *                          Structure: Sparse (CSC) Matrix (numRows * numCols).
 *                          Type: 'Double'.
 *                          Range: (-inf, inf).
 * References:
 *	1. 	https://github.com/pymatting/pymatting/blob/master/pymatting/preconditioner/ichol.py
 *	2.	Incomplete Cholesky Decomposition: https://en.wikipedia.org/wiki/Incomplete_Cholesky_factorization
 * Remarks:
 *	1.	The Sparse Matrices are given in Compressed Sparse Column (CSC) format.
 * TODO:
 *	1.	Add "Zero Fill" variant of the decomposition.
 * Release Notes:
 *	-	1.0.000	10/07/2020	Royi Avital
 *		*	First release version.
 */


#include <stdlib.h>

#include "mex.h"
#include "IncompleteCholeskyDecomposition.c"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	mwSize numRows, numCols, maxNumNzA, numNzA, numNz;
	mwIndex *Ir, *Jc; // Pseudo Row / Column Vectors of Sparse Matrix
	double *vData, *vDataA, discardThr, shiftVal;
	unsigned int *vIndices, *vIndicesPtr, *vIndicesA, *vIndicesPtrA, numShifts, maxNumNz, ii;
	int icholNnz;
	
	// Validate Input
	if ( nrhs != 4 )
	{
		mexErrMsgIdAndTxt("IncompleteCholeskyDecomposition:IncompleteCholeskyDecomposition:nrhs", "There must be 4 inputs: A real sparse positive definite matrix of type double, and 3 scalars of type double");
	}
	// Validating the 1st input
	if( !mxIsNumeric(prhs[0]) || mxIsEmpty(prhs[0]) || mxIsComplex(prhs[0]) || (mxGetClassID(prhs[0]) != mxDOUBLE_CLASS) || !mxIsSparse(prhs[0]) || (mxGetM(prhs[0]) != mxGetN(prhs[0])) )
	{
		mexErrMsgIdAndTxt("IncompleteCholeskyDecomposition:IncompleteCholeskyDecomposition:nrhs", "The 1st input must be a real sparse positive definite matrix of type double");
	}

	// Validating the 2nd input
	if ( !mxIsNumeric(prhs[1]) || mxIsEmpty(prhs[1]) || mxIsComplex(prhs[1]) || (mxGetClassID(prhs[1]) != mxDOUBLE_CLASS) || (mxGetNumberOfElements(prhs[1]) != 1) )
	{
		mexErrMsgIdAndTxt("IncompleteCholeskyDecomposition:IncompleteCholeskyDecomposition:nrhs", "Discarding Threshold: The 2nd input must be a real non negative scalar of type double");
	}

	// Validating the 3rd input
	if ( !mxIsNumeric(prhs[2]) || mxIsEmpty(prhs[2]) || mxIsComplex(prhs[2]) || (mxGetClassID(prhs[2]) != mxDOUBLE_CLASS) || (mxGetNumberOfElements(prhs[2]) != 1) )
	{
		mexErrMsgIdAndTxt("IncompleteCholeskyDecomposition:IncompleteCholeskyDecomposition:nrhs", "Shift Value: The 3rd input must be a real non negative scalar of type double");
	}

	// Validating the 4th input
	if ( !mxIsNumeric(prhs[3]) || mxIsEmpty(prhs[3]) || mxIsComplex(prhs[3]) || (mxGetClassID(prhs[3]) != mxDOUBLE_CLASS) || (mxGetNumberOfElements(prhs[3]) != 1) )
	{
		mexErrMsgIdAndTxt("IncompleteCholeskyDecomposition:IncompleteCholeskyDecomposition:nrhs", "Maximum Number of Non Zeros: The 4th input must be a real positive scalar of type double");
	}

	// Validate Output
	if (nlhs != 1)
	{
		mexErrMsgIdAndTxt("IncompleteCholeskyDecomposition:IncompleteCholeskyDecomposition:nlhs", "Output must be double real sparse matrix");
	}
	
	// Input Matrix Dimensions
	numRows = mxGetM(prhs[0]);
	numCols = mxGetN(prhs[0]);
	
	// Input Matrix Data
	vDataA = (double*)mxGetData(prhs[0]);
	// Allocated on creation of the Sparse Matrix, Length as the number of non zero elements
	Ir = mxGetIr(prhs[0]);
	// Allocated on creation of the Sparse Matrix, Length as the number of columns + 1
	Jc = mxGetJc(prhs[0]);

	numNzA = Jc[numCols]; // The last element of Jc has the actual number of NNZ in the Matrix.
	maxNumNzA = mxGetNzmax(prhs[0]); // The maximum number of elements in mA (Length of Ir and vDataA).

	discardThr	= (double)mxGetScalar(prhs[1]);
	shiftVal	= (double)mxGetScalar(prhs[2]);
	maxNumNz	= (unsigned int)mxGetScalar(prhs[3]);

	if (discardThr < 0.0)
	{
		mexErrMsgIdAndTxt("IncompleteCholeskyDecomposition:IncompleteCholeskyDecomposition:nrhs", "Discarding Threshold: The 2nd input must be a real non negative scalar of type double");
	}

	if (shiftVal < 0.0)
	{
		mexErrMsgIdAndTxt("IncompleteCholeskyDecomposition:IncompleteCholeskyDecomposition:nrhs", "Shift Value: The 3rd input must be a real non negative scalar of type double");
	}

	vData		= (double *)malloc(maxNumNz * sizeof(double));
	vIndices	= (unsigned int *)malloc(maxNumNz * sizeof(unsigned int));
	vIndicesPtr = (unsigned int*)calloc(numCols + 1, sizeof(unsigned int));

	vIndicesA		= (unsigned int*)malloc(numNzA * sizeof(unsigned int));
	vIndicesPtrA	= (unsigned int*)malloc((numCols + 1) * sizeof(unsigned int));

	for (ii = 0; ii < numNzA; ii++)
	{
		vIndicesA[ii] = (unsigned int)Ir[ii];
	}

	for (ii = 0; ii < numCols + 1; ii++)
	{
		vIndicesPtrA[ii] = (unsigned int)Jc[ii];
	}
	
	// IncompleteCholeskyDecomposition(vData, vIndices, vIndicesPtr, vDataA, vIndicesA, vIndicesPtrA, (unsigned int)numCols, discardThr, vShifts, numShifts, maxNumNz);
	icholNnz = _IncompleteCholDec(vData, vIndices, vIndicesPtr, vDataA, vIndicesA, vIndicesPtrA, (unsigned int)numCols, discardThr, shiftVal, maxNumNz);

	if (icholNnz == -1)
	{
		free(vData);
		free(vIndices);
		free(vIndicesPtr);

		free(vIndicesA);
		free(vIndicesPtrA);
		mexErrMsgIdAndTxt("IncompleteCholeskyDecomposition:IncompleteCholeskyDecomposition:nlhs", "Thresholded incomplete Cholesky decomposition failed due to insufficient positive definiteness of the matrix A.\n");
	}

	if (icholNnz == -2)
	{
		free(vData);
		free(vIndices);
		free(vIndicesPtr);

		free(vIndicesA);
		free(vIndicesPtrA);
		mexErrMsgIdAndTxt("IncompleteCholeskyDecomposition:IncompleteCholeskyDecomposition:nlhs", "Thresholded incomplete Cholesky decomposition failed because more than maxNumNz non zero elements were created.\n");
	}

	numNz = (mwSize)vIndicesPtr[numCols];

	//// Allocating the output Sparse Matrix
	plhs[0] = mxCreateSparse(numRows, numCols, numNz, mxREAL);

	vDataA = (double*)mxGetData(plhs[0]);

	// Allocated on creation of the Sparse Matrix, Length as the number of non zero elements
	Ir = mxGetIr(plhs[0]);
	// Allocated on creation of the Sparse Matrix, Length as the number of columns + 1
	Jc = mxGetJc(plhs[0]);

	for (ii = 0; ii < numNz; ii++)
	{
		vDataA[ii] = vData[ii];
		//mexPrintf("vData[%u] = %f\n", ii, vData[ii]);
		//mexPrintf("vDataA[%u] = %f\n", ii, vDataA[ii]);
		Ir[ii] = vIndices[ii];
		//mexPrintf("vIndices[%u] = %u\n", ii, vIndices[ii]);
		//mexPrintf("Ir[%u] = %u\n", ii, (unsigned int)Ir[ii]);
	}

	for (ii = 0; ii < numCols + 1; ii++)
	{
		Jc[ii] = vIndicesPtr[ii];
		//mexPrintf("vIndicesPtr[%u] = %u\n", ii, vIndicesPtr[ii]);
		//mexPrintf("Jc[%u] = %u\n", ii, (unsigned int)Jc[ii]);
	}

	free(vData);
	free(vIndices);
	free(vIndicesPtr);

	free(vIndicesA);
	free(vIndicesPtrA);
	
	
}