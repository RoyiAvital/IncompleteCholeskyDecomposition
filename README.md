# Incomplete Cholesky Decomposition Threshold

Implementation of the Incomplete Cholesky Decomposition with Thresholding.  
The project includes a `C` implementation with MATLAB MEX wrapper.

## Generating MATLAB MEX

 1.	Download the repository.
 2. Run `MakeMex` in MATLAB with pre defined MATLAB MEX Compiler.
 3.	Go through the Unit Tests and the Run Time Analysis.
 
The MEX Wrapper supports only Sparse Real Matrices of Type Double.

## Performance

Comparing the performance with MATLAB's functions.

## Decomposition

![](https://i.imgur.com/zYKNq9o.png)

The MEX file and MATLAB's `ICT` were the most memory efficient.

## Pre Conditioning (Solving the Linear System)
 
![](https://i.imgur.com/iwoVKdM.png)

## References

 * 	[PyMatting](https://github.com/pymatting/pymatting).  
	The `C` code is basically a redo of the Pre Conditioner in the Python package.
 *	[MATLAB `ichol()`](https://www.mathworks.com/help/matlab/ref/ichol.html).