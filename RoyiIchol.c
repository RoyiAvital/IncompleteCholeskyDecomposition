#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <assert.h>

// Code by Thomas Germer
// https://github.com/99991/testing
// Used in order to compare against MATLAB after my results showed my code is slower

#ifdef RUN_TIME
double sec(){
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + 1e-9 * t.tv_nsec;
}
#endif

#define NEW(n, T) calloc((n), sizeof(T))
#define FOR(i, a, b) for (int i = (a); i < (b); i++)

void _backsub_L_csc_inplace(const double *data, const int *indices, const int *indptr, double *x, int n){
    for (int j = 0; j < n; j++){
        int k = indptr[j];
        double L_jj = data[k];
        double temp = x[j] / L_jj;

        x[j] = temp;

        for (int k = indptr[j] + 1; k < indptr[j + 1]; k++){
            int i = indices[k];
            double L_ij = data[k];

            x[i] -= L_ij * temp;
        }
    }
}

void _backsub_LT_csc_inplace(const double *data, const int *indices, const int *indptr, double *x, int n){
    for (int i = n - 1; i >= 0; i--){
        double s = x[i];

        for (int k = indptr[i] + 1; k < indptr[i + 1]; k++){
            int j = indices[k];
            double L_ji = data[k];
            s -= L_ji * x[j];
        }
        int k = indptr[i];
        double L_ii = data[k];

        x[i] = s / L_ii;
    }
}

void mergesort(int *x, size_t n, int *tmp){
    if (n <= 10){
        for (size_t i = 1; i < n; i++){
            int value = x[i];
            size_t j;
            for (j = i; j > 0 && x[j - 1] > value; j--) x[j] = x[j - 1];
            x[j] = value;
        }
    }else{
        size_t m = n / 2, i = 0, j = m, c = 0;
        mergesort(x, m, tmp);
        mergesort(x + m, n - m, tmp);
        while (i < m && j < n) tmp[c++] = x[i] < x[j] ? x[i++] : x[j++];
        while (i < m) tmp[c++] = x[i++];
        for (size_t k = 0; k < c; k++) x[k] = tmp[k];
    }
}

int _ichol(
    int n,
    const double *Av,
    const int *Ar,
    const int *Ap,
    double *Lv,
    int *Lr,
    int *Lp,
    double discard_threshold,
    double relative_threshold,
    double shift,
    int max_nnz
){
    int nnz = 0;
    int c_n = 0;
    int *s = NEW(n, int);  // Next non-zero row index i in column j of L
    int *t = NEW(n, int);  // First subdiagonal index i in column j of A
    int *l = NEW(n, int);  // Linked list of non-zero columns in row k of L
    double  *a = NEW(n, double );  // Values of column j
    uint8_t *b = NEW(n, uint8_t);  // b[i] indicates if the i-th element of column j is non-zero
    int *c = NEW(n, int);  // Row indices of non-zero elements in column j
    double  *d = NEW(n, double );  // Diagonal elements of A
    int *e = NEW(n, int);  // Temporary array for mergesort
    FOR(i, 0, n){
        l[i] = -1;
        d[i] = shift;
    }

    FOR(j, 0, n){
        FOR(idx, Ap[j], Ap[j + 1]){
            int i = Ar[idx];
            if (i == j){
                d[j] += Av[idx];
                t[j] = idx + 1;
            }
        }
    }
	
#ifdef RUN_TIME
    double dt[4] = {0.0, 0.0, 0.0, 0.0};
#endif

    FOR(j, 0, n){  // For each column j

#ifdef RUN_TIME
        double t0 = sec();
#endif

        double l1_norm_column = 0.0;
        FOR(idx, t[j], Ap[j + 1]){  // For each L_ij
            int i = Ar[idx];
            double L_ij = Av[idx];
            l1_norm_column += fabs(L_ij);

            if (L_ij != 0.0 && i > j){
                a[i] += L_ij;  // Assign non-zero value to L_ij in sparse column
                if (!b[i]){
                    b[i] = 1;  // Mark it as non-zero
                    c[c_n] = i;  // Remember index for later deletion
                    c_n += 1;
                }
            }
        }

#ifdef RUN_TIME
        double t1 = sec();
#endif
        int k = l[j];  // Find index k of column with non-zero element in row j
        while (k != -1){  // For each column of that type
            int k0 = s[k];  // Start index of non-zero elements in column k
            int k1 = Lp[k + 1];  // End index
            int k2 = l[k];  // Remember next column index before it is overwritten
            double L_jk = Lv[k0];  // Value of non-zero element at start of column
            k0 += 1;  // Advance to next non-zero element in column
            if (k0 < k1){  // If there is a next non-zero element
                s[k] = k0;  // Advance start index in column k to next non-zero element
                int i = Lr[k0];  // Row index of next non-zero element in column k
                l[k] = l[i];  // Remember old list i index in list k
                l[i] = k;  // Insert index of non-zero element into list i
                FOR(idx, k0, k1){  // For each non-zero L_ik in column k
                    int i = Lr[idx];
                    double L_ik = Lv[idx];
                    a[i] -= L_ik * L_jk;  // Update element L_ij in sparse column
                    if (!b[i]){  // Check if sparse column element was zero
                        b[i] = 1;  // Mark as non-zero in sparse column
                        c[c_n] = i;  // Remember index for later deletion
                        c_n += 1;
                    }
                }
            }
            k = k2;  // Advance to next column k
        }

        if (d[j] <= 0.0) return -1;
        if (nnz + 1 + c_n > max_nnz) return -2;
        d[j] = sqrt(d[j]);  // Update diagonal element L_ii
        Lv[nnz] = d[j];  // Add diagonal element L_ii to L
        Lr[nnz] = j;  // Add row index of L_ii to L
        nnz += 1;
        s[j] = nnz;  // Set first non-zero index of column j
#ifdef RUN_TIME
        double t2 = sec();
#endif
        mergesort(c, c_n, e);  // Sort row indices of column j for correct insertion order into L

#ifdef RUN_TIME
        double t3 = sec();
#endif
        FOR(idx, 0, c_n){
            int i = c[idx];
            double L_ij = a[i] / d[j];  // Get non-zero element from sparse column j
            d[i] -= L_ij * L_ij;  // Update diagonal element L_ii
            if (i == j || ((fabs(L_ij) > discard_threshold) && (fabs(L_ij) > (l1_norm_column * relative_threshold)))){  // If element is sufficiently non-zero or diagonal
                Lv[nnz] = L_ij;  // Add element L_ij to L
                Lr[nnz] = i;  // Add row index of L_ij
                nnz += 1;
            }
            a[i] = 0.0;  // Set element i in column j to zero
            b[i] = 0;  // Mark element as zero
        }
        c_n = 0;  // Discard row indices of non-zero elements in column j.
        Lp[j + 1] = nnz;  // Update count of non-zero elements up to column j
        if (Lp[j] + 1 < Lp[j + 1]){  // If column j has a non-zero element below diagonal
            int i = Lr[Lp[j] + 1];  // Row index of first off-diagonal non-zero element
            l[j] = l[i];  // Remember old list i index in list j
            l[i] = j;  // Insert index of non-zero element into list i
        }
#ifdef RUN_TIME
        double t4 = sec();

        dt[0] += t1 - t0;
        dt[1] += t2 - t1;
        dt[2] += t3 - t2;
        dt[3] += t4 - t3;
#endif
    }
#ifdef RUN_TIME
    printf("reading    : %f sec\n", dt[0]);
    printf("multiplying: %f sec\n", dt[1]);
    printf("sorting    : %f sec\n", dt[2]);
    printf("writing    : %f sec\n", dt[3]);
    printf("total      : %f sec\n", dt[0] + dt[1] + dt[2] + dt[3]);
#endif
    free(s);
    free(t);
    free(l);
    free(a);
    free(b);
    free(c);
    free(d);
    free(e);

    return nnz;
}
