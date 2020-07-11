function [ mA ] = GenWlsMatrix( numRows )
%GenWlsMatrix generates the Weighted LEast Squares operator in Sparse
%Matrix form.
%   Detailed explanation goes here

paramAlpha  = 1.2;
paramLambda = 0.5;

mI = rand(numRows, numRows);

EPS_NUM     = eps();
SMALL_NUM   = 1e-4;

numRows     = size(mI, 1);
numCols     = size(mI, 2);
numPixel    = numRows * numCols;

mL = log(mI + EPS_NUM);

% Compute affinities between adjacent pixels based on gradients of L
% mDy = mL(2:end, 2:end) - mL(1:end - 1, 2:end);
mDy = diff(mL, 1, 1);
mDy = -paramLambda ./ ((abs(mDy) .^ paramAlpha) + SMALL_NUM);
mDy = padarray(mDy, [1, 0], 'post');

% mDx = mL(2:end, 2:end) - mL(2:end - 1, 1:end - 1);
mDx = diff(mL, 1, 2); 
mDx = -paramLambda ./ ((abs(mDx) .^ paramAlpha) + SMALL_NUM);
mDx = padarray(mDx, [0, 1], 'post');

% Construct a five point spatially inhomogeneous Laplacian matrix
mB = [mDx(:), mDy(:)];
vDiagIdx = [-numRows, -1];
mA = spdiags(mB, vDiagIdx, numPixel, numPixel);

mA = mA + mA.';
mA = mA + spdiags(1 - sum(mA, 2), 0, numPixel, numPixel);


end

