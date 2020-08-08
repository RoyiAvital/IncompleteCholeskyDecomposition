
rng(0);

numRows     = 300;
discardThr  = 1e-4;
maxNumNz    = 5e6;

sIchol = struct('type', 'ict', 'droptol', discardThr / 1000);

mA = sprandn(numRows, numRows, 0.9);
mA = (mA.' * mA) + speye(numRows);
mB = full(mA);

vX = ones(numRows, 1);

mL = ichol(mA, sIchol);
disp(['Error Norm - ', num2str(norm(mA * vX - (mL * mL.') * vX, 'fro'))]);
mLA = full(mL);
mL = IncompleteCholeskyDecompositionMex(mA, discardThr, maxNumNz);
disp(['Error Norm - ', num2str(norm(mA * vX - (mL * mL.') * vX, 'fro'))]);
mLB = full(mL);

vX = randn(numRows, 1);

vY = SolveCholsekyLinearSystemMex(mL, vX);
vYY = (mL * mL.') \ vX;
% SolveCholsekyLinearSystemMex(mL, vX);
% vYYY = mA \ vX;

norm(vY - vYY, 'fro')
