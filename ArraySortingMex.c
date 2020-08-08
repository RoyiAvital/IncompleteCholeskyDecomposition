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

#define ELEMENT_TYPE_IDX 1

#include <stdlib.h>

#include "mex.h"
#include "ArraySorting.c"

// See https://stackoverflow.com/questions/2335888

#if ELEMENT_TYPE_IDX == 1
#define MX_TYPE mxUINT32_CLASS
#elif ELEMENT_TYPE_IDX == 2
#define MX_TYPE mxINT32_CLASS
#elif ELEMENT_TYPE_IDX == 3
#define MX_TYPE mxSINGLE_CLASS
#elif ELEMENT_TYPE_IDX == 4
#define MX_TYPE mxDOUBLE_CLASS
#endif // ELEMENT_TYPE_IDX == 1


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	mwSize numElements;
#if ELEMENT_TYPE_IDX == 1
	unsigned int *vA;
#elif ELEMENT_TYPE_IDX == 2
	int *vA;
#elif ELEMENT_TYPE_IDX == 3
	single *vA;
#elif ELEMENT_TYPE_IDX == 4
	double *vA;
#endif // ELEMENT_TYPE_IDX == 1

	
	// Validating Input
	if( (nrhs != 1) || !mxIsNumeric(prhs[0]) || mxIsEmpty(prhs[0]) || mxIsComplex(prhs[0]) || (mxGetClassID(prhs[0]) != MX_TYPE) )
	{
#if ELEMENT_TYPE_IDX == 1
		mexErrMsgIdAndTxt("ArraySorting:ArraySorting:nrhs", "The input must be a single real numeric array input of type UINT32");
#elif ELEMENT_TYPE_IDX == 2
		mexErrMsgIdAndTxt("ArraySorting:ArraySorting:nrhs", "The input must be a single real numeric array input of type INT32");
#elif ELEMENT_TYPE_IDX == 3
		mexErrMsgIdAndTxt("ArraySorting:ArraySorting:nrhs", "The input must be a single real numeric array input of type Single");
#elif ELEMENT_TYPE_IDX == 4
		mexErrMsgIdAndTxt("ArraySorting:ArraySorting:nrhs", "The input must be a single real numeric array input of type Double");
#endif // ELEMENT_TYPE_IDX == 1
	}

	// Validate Output
	if (nlhs != 0)
	{
		mexErrMsgIdAndTxt("ArraySorting:ArraySorting:nlhs", "There must be no output. The function is in place.");
	}
	
	// Input Matrix Dimensions
	numElements = mxGetNumberOfElements(prhs[0]);
	
	// Input Matrix Data
#if ELEMENT_TYPE_IDX == 1
	vA = (unsigned int*)mxGetData(prhs[0]);
#elif ELEMENT_TYPE_IDX == 2
	vA = (int*)mxGetData(prhs[0]);
#elif ELEMENT_TYPE_IDX == 3
	vA = (single*)mxGetData(prhs[0]);
#elif ELEMENT_TYPE_IDX == 4
	vA = (double*)mxGetData(prhs[0]);
#endif // ELEMENT_TYPE_IDX == 1
	
	ArrayQuickSort(vA, (unsigned int)0, (unsigned int)numElements - (unsigned int)1);
	
	
}