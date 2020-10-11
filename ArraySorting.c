/*
 * Array Sorting
 * Sorts a contiguous array using Sorting Network and High Level sorting algorithm.
 *
 * References:
 *	1. 	Based on https://github.com/eggpi/sorting-networks-test by Guilherme P. Goncalves.
 * Remarks:
 *	1.	Supports `unisgned int`, `int`, `single` and `double` by changing `ELEMENT_TYPE_IDX`.
 * TODO:
 *	1.	Add option for Merge Sort in addition to Quick Sort.
 * Release Notes:
 *	-	1.0.000	13/07/2020	Royi Avital
 *		*	First release version.
 */

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

#define SWAP(valA, valB) \
	{ELEMENT_TYPE tmpVal; \
	(tmpVal) = (valA); \
	(valA) = (valB); \
	(valB) = (tmpVal); \
}

#define SWAP_PTR(valA, valB) \
	{ELEMENT_TYPE *tmpVal; \
	(tmpVal) = (valA); \
	(valA) = (valB); \
	(valB) = (tmpVal); \
}


#define SWAP_SORT(vA, ii, jj) \
if ((vA)[ii] > (vA)[jj]) { \
	ELEMENT_TYPE tmpVal; \
	tmpVal = (vA)[ii]; \
	(vA)[ii] = (vA)[jj]; \
	(vA)[jj] = tmpVal; \
}

#define NSORT2(vA) \
SWAP_SORT((vA), 0, 1)

#define NSORT3(vA) \
SWAP_SORT((vA), 1, 2); SWAP_SORT((vA), 0, 2); SWAP_SORT((vA), 0, 1);

#define NSORT4(vA) \
SWAP_SORT((vA), 0, 1); SWAP_SORT((vA), 2, 3); SWAP_SORT((vA), 0, 2); \
SWAP_SORT((vA), 1, 3); SWAP_SORT((vA), 1, 2);

#define NSORT5(vA) \
SWAP_SORT((vA), 0, 1); SWAP_SORT((vA), 3, 4); SWAP_SORT((vA), 2, 4); \
SWAP_SORT((vA), 2, 3); SWAP_SORT((vA), 0, 3); SWAP_SORT((vA), 0, 2); \
SWAP_SORT((vA), 1, 4); SWAP_SORT((vA), 1, 3); SWAP_SORT((vA), 1, 2);

#define NSORT6(vA) \
SWAP_SORT((vA), 1, 2); SWAP_SORT((vA), 0, 2); SWAP_SORT((vA), 0, 1); \
SWAP_SORT((vA), 4, 5); SWAP_SORT((vA), 3, 5); SWAP_SORT((vA), 3, 4); \
SWAP_SORT((vA), 0, 3); SWAP_SORT((vA), 1, 4); SWAP_SORT((vA), 2, 5); \
SWAP_SORT((vA), 2, 4); SWAP_SORT((vA), 1, 3); SWAP_SORT((vA), 2, 3);

#define NSORT7(vA) \
SWAP_SORT((vA), 1, 2); SWAP_SORT((vA), 0, 2); SWAP_SORT((vA), 0, 1); \
SWAP_SORT((vA), 3, 4); SWAP_SORT((vA), 5, 6); SWAP_SORT((vA), 3, 5); \
SWAP_SORT((vA), 4, 6); SWAP_SORT((vA), 4, 5); SWAP_SORT((vA), 0, 4); \
SWAP_SORT((vA), 0, 3); SWAP_SORT((vA), 1, 5); SWAP_SORT((vA), 2, 6); \
SWAP_SORT((vA), 2, 5); SWAP_SORT((vA), 1, 3); SWAP_SORT((vA), 2, 4); \
SWAP_SORT((vA), 2, 3);

#define NSORT8(vA) \
SWAP_SORT((vA), 0, 1); SWAP_SORT((vA), 2, 3); SWAP_SORT((vA), 0, 2); \
SWAP_SORT((vA), 1, 3); SWAP_SORT((vA), 1, 2); SWAP_SORT((vA), 4, 5); \
SWAP_SORT((vA), 6, 7); SWAP_SORT((vA), 4, 6); SWAP_SORT((vA), 5, 7); \
SWAP_SORT((vA), 5, 6); SWAP_SORT((vA), 0, 4); SWAP_SORT((vA), 1, 5); \
SWAP_SORT((vA), 1, 4); SWAP_SORT((vA), 2, 6); SWAP_SORT((vA), 3, 7); \
SWAP_SORT((vA), 3, 6); SWAP_SORT((vA), 2, 4); SWAP_SORT((vA), 3, 5); \
SWAP_SORT((vA), 3, 4);


int ArrayPartition(ELEMENT_TYPE *vA, unsigned int s, unsigned int e) {
	unsigned int ii, jj, p;
	ELEMENT_TYPE tmpVal;

	/* Get pivot from median of three */
	p = (s + e)/2;

	if (vA[p] > vA[s]) {
		if (vA[e] > vA[s]) {
			/* p > e > s or e > p > s */
			p = (vA[p] > vA[e]) ? e : p;
		} else p = s;
	} else {
		if (vA[e] > vA[p]) {
			/* s > e > p or e > s > p */
			p = (vA[e] > vA[s]) ? s : e;
		}
	}

	/* Place pivot in the beginning */
	tmpVal = vA[p];
	vA[p] = vA[s];
	vA[s] = tmpVal;

	/* Partition <= | pivot | >= */
	ii = s;
	jj = e;
	while (1) {
		while (ii < e && vA[ii] <= vA[s]) ii++;
		while (jj > s && vA[jj] >= vA[s]) jj--;

		if (ii > jj) break;

		tmpVal = vA[ii];
		vA[ii] = vA[jj];
		vA[jj] = tmpVal;

		ii++;
		jj--;
	}

	/* Swap pivot back to its place */
	tmpVal = vA[--ii];
	vA[ii] = vA[s];
	vA[s] = tmpVal;

	return ii;
}

static void ArrayQuickSort(ELEMENT_TYPE *vA, unsigned int s, unsigned int e) {
	unsigned int p, size;

	size = e - s + 1;
	if (size <= 1) {
		return;
	}

	switch (size) {
		case 2:
			NSORT2(vA + s);
			break;
		case 3:
			NSORT3(vA + s);
			break;
		case 4:
			NSORT4(vA + s);
			break;
		case 5:
			NSORT5(vA + s);
			break;
		case 6:
			NSORT6(vA + s);
			break;
		case 7:
			NSORT7(vA + s);
			break;
		case 8:
			NSORT8(vA + s);
			break;
		default:
			p = ArrayPartition(vA, s, e);
			ArrayQuickSort(vA, s, p - 1);
			ArrayQuickSort(vA, p + 1, e);
	}

	return;
}

void ArrayMergeSort(ELEMENT_TYPE *x, unsigned int numElements, ELEMENT_TYPE *vTmp) {
	if (numElements <= 10) {
		for (unsigned int ii = 1; ii < numElements; ii++) {
			ELEMENT_TYPE tmpVal = x[ii];
			unsigned int jj;
			for (jj = ii; jj > 0 && x[jj - 1] > tmpVal; jj--) x[jj] = x[jj - 1];
			x[jj] = tmpVal;
		}
	}
	else {
		unsigned int m = numElements / 2, ii = 0, jj = m, ll = 0;
		ArrayMergeSort(x, m, vTmp);
		ArrayMergeSort(x + m, numElements - m, vTmp);
		while (ii < m && jj < numElements) vTmp[ll++] = x[ii] < x[jj] ? x[ii++] : x[jj++];
		while (ii < m) vTmp[ll++] = x[ii++];
		for (unsigned int kk = 0; kk < ll; kk++) x[kk] = vTmp[kk];
	}
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
        SWAP_PTR(a, b); // swap ptrs
        rsz <<= 2;     // quadruple run size
    }
    return;                           // return sorted array
}