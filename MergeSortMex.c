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

#define ELEMENT_TYPE_IDX 4

#include <stdlib.h>
#include "mex.h"

#ifndef ELEMENT_TYPE_IDX
#define ELEMENT_TYPE_IDX 1
#endif // !ELEMENT_TYPE_IDX


#if ELEMENT_TYPE_IDX == 1
#define ELEMENT_TYPE unsigned int
#elif ELEMENT_TYPE_IDX == 2
#define ELEMENT_TYPE int
#elif ELEMENT_TYPE_IDX == 3
#define ELEMENT_TYPE float
#elif ELEMENT_TYPE_IDX == 4
#define ELEMENT_TYPE double
#endif // ELEMENT_TYPE_IDX == 1
 
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

#define SWAP(valA, valB) \
	{ELEMENT_TYPE tmpVal; \
	(tmpVal) = (valA); \
	(valA) = (valB); \
	(valB) = (tmpVal); \
}

#define SWAPPTR(valA, valB) \
	{ELEMENT_TYPE *tmpVal; \
	(tmpVal) = (valA); \
	(valA) = (valB); \
	(valB) = (tmpVal); \
}


void BottomUpMergeSort(ELEMENT_TYPE a[], ELEMENT_TYPE b[], unsigned int n)
{
    ELEMENT_TYPE* p0r;       // ptr to      run 0
    ELEMENT_TYPE* p0e;       // ptr to end  run 0
    ELEMENT_TYPE* p1r;       // ptr to      run 1
    ELEMENT_TYPE* p1e;       // ptr to end  run 1
    ELEMENT_TYPE* p2r;       // ptr to      run 2
    ELEMENT_TYPE* p2e;       // ptr to end  run 2
    ELEMENT_TYPE* p3r;       // ptr to      run 3
    ELEMENT_TYPE* p3e;       // ptr to end  run 3
    ELEMENT_TYPE* pax;       // ptr to set of runs in a
    ELEMENT_TYPE* pbx;       // ptr for merged output to b
    unsigned int rsz = 1; // run size
    if (n < 2)
        return;
    if (n == 2) {
        // if (a[0] > a[1])std::swap(a[0], a[1]);
        if (a[0] > a[1]) { SWAP(a[0], a[1]) };
        return;
    }
    if (n == 3) {
        // if (a[0] > a[2])std::swap(a[0], a[2]);
        // if (a[0] > a[1])std::swap(a[0], a[1]);
        // if (a[1] > a[2])std::swap(a[1], a[2]);
        if (a[0] > a[2]) { SWAP(a[0], a[2]) };
        if (a[0] > a[1]) { SWAP(a[0], a[1]) };
        if (a[1] > a[2]) { SWAP(a[1], a[2]) };
        return;
    }
    while (rsz < n) {
        pbx = &b[0];
        pax = &a[0];
        while (pax < &a[n]) {
            p0e = rsz + (p0r = pax);
            if (p0e >= &a[n]) {
                p0e = &a[n];
                goto cpy10;
            }
            p1e = rsz + (p1r = p0e);
            if (p1e >= &a[n]) {
                p1e = &a[n];
                goto mrg201;
            }
            p2e = rsz + (p2r = p1e);
            if (p2e >= &a[n]) {
                p2e = &a[n];
                goto mrg3012;
            }
            p3e = rsz + (p3r = p2e);
            if (p3e >= &a[n])
                p3e = &a[n];
            // 4 way merge
            while (1) {
                if (*p0r <= *p1r) {
                    if (*p2r <= *p3r) {
                        if (*p0r <= *p2r) {
                        mrg40:                      *pbx++ = *p0r++;    // run 0 smallest
                            if (p0r < p0e)       // if not end run continue
                                continue;
                            goto mrg3123;       // merge 1,2,3
                        }
                        else {
                        mrg42:                      *pbx++ = *p2r++;    // run 2 smallest
                            if (p2r < p2e)       // if not end run continue
                                continue;
                            goto mrg3013;       // merge 0,1,3
                        }
                    }
                    else {
                        if (*p0r <= *p3r) {
                            goto mrg40;         // run 0 smallext
                        }
                        else {
                        mrg43:                      *pbx++ = *p3r++;    // run 3 smallest
                            if (p3r < p3e)       // if not end run continue
                                continue;
                            goto mrg3012;       // merge 0,1,2
                        }
                    }
                }
                else {
                    if (*p2r <= *p3r) {
                        if (*p1r <= *p2r) {
                        mrg41:                      *pbx++ = *p1r++;    // run 1 smallest
                            if (p1r < p1e)       // if not end run continue
                                continue;
                            goto mrg3023;       // merge 0,2,3
                        }
                        else {
                            goto mrg42;         // run 2 smallest
                        }
                    }
                    else {
                        if (*p1r <= *p3r) {
                            goto mrg41;         // run 1 smallest
                        }
                        else {
                            goto mrg43;         // run 3 smallest
                        }
                    }
                }
            }
            // 3 way merge
        mrg3123:    p0r = p1r;
            p0e = p1e;
        mrg3023:    p1r = p2r;
            p1e = p2e;
        mrg3013:    p2r = p3r;
            p2e = p3e;
        mrg3012:    while (1) {
            if (*p0r <= *p1r) {
                if (*p0r <= *p2r) {
                    *pbx++ = *p0r++;        // run 0 smallest
                    if (p0r < p0e)           // if not end run continue
                        continue;
                    goto mrg212;            // merge 1,2
                }
                else {
                mrg32:                  *pbx++ = *p2r++;        // run 2 smallest
                    if (p2r < p2e)           // if not end run continue
                        continue;
                    goto mrg201;            // merge 0,1
                }
            }
            else {
                if (*p1r <= *p2r) {
                    *pbx++ = *p1r++;        // run 1 smallest
                    if (p1r < p1e)           // if not end run continue
                        continue;
                    goto mrg202;            // merge 0,2
                }
                else {
                    goto mrg32;             // run 2 smallest
                }
            }
        }
        // 2 way merge
    mrg212:     p0r = p1r;
        p0e = p1e;
    mrg202:     p1r = p2r;
        p1e = p2e;
    mrg201:     while (1) {
        if (*p0r <= *p1r) {
            *pbx++ = *p0r++;            // run 0 smallest
            if (p0r < p0e)               // if not end run continue
                continue;
            goto cpy11;
        }
        else {
            *pbx++ = *p1r++;            // run 1 smallest
            if (p1r < p1e)               // if not end run continue
                continue;
            goto cpy10;
        }
    }
    // 1 way copy
cpy11:      p0r = p1r;
    p0e = p1e;
cpy10:      while (1) {
    *pbx++ = *p0r++;                // copy element
    if (p0r < p0e)                  // if not end of run continue
        continue;
    break;
}
pax += rsz << 2;            // setup for next set of runs
        }
        // std::swap(a, b);                // swap ptrs
        SWAPPTR(a, b); // swap ptrs
        rsz <<= 2;     // quadruple run size
    }
    return;                           // return sorted array
}


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	mwSize numElements;
#if ELEMENT_TYPE_IDX == 1
	unsigned int *vA, *vB;
#elif ELEMENT_TYPE_IDX == 2
	int *vA, *vB;
#elif ELEMENT_TYPE_IDX == 3
	single *vA, *vB;
#elif ELEMENT_TYPE_IDX == 4
	double *vA, *vB;
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
    vA = (double*)mxGetData(prhs[0]);
    vB = (ELEMENT_TYPE*)malloc(numElements * sizeof(ELEMENT_TYPE));
	
    BottomUpMergeSort(vA, vB, numElements);
	
}