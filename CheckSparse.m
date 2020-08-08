

numRows = 3;
numCols = 3;

%{
From SciPy:
import numpy as np
from scipy.sparse import csc_matrix

mA = np.array([[1, 0, 0], [0, 2, 0], [0, 0, 3]])

mA
Out[4]: 
array([[1, 0, 0],
       [0, 2, 0],
       [0, 0, 3]])

mB = csc_matrix(mA);

print(mB)
  (0, 0)	1
  (1, 1)	2
  (2, 2)	3

mB.indices
Out[7]: array([0, 1, 2], dtype=int32)

mB.indptr
Out[8]: array([0, 1, 2, 3], dtype=int32)
%}
vV = [1; 2; 3];
vR = [0; 1; 2];
vP = [0; 1; 2; 3];

mA = BuildCscSparseMatrix(numRows, numCols, vV, vR, vP);
