
rng(0);

numRows     = 1000;
discardThr  = 1e-6;
relThr      = 1e-4;
shiftVal    = 0;
maxNumNz    = 5e6;

sIchol = struct('type', 'ict', 'droptol', relThr);

mA = sprandn(numRows, numRows, 0.15);
mA = (mA.' * mA) + speye(numRows);
mB = full(mA);

vX = ones(numRows, 1);

mL = ichol(mA, sIchol);
disp(['Error Norm - ', num2str(norm(mA * vX - (mL * mL.') * vX, 'fro'))]);
mLA = full(mL);
mL = IncompleteCholeskyDecompositionMex(mA, discardThr, shiftVal, maxNumNz);
disp(['Error Norm - ', num2str(norm(mA * vX - (mL * mL.') * vX, 'fro'))]);
mLB = full(mL);
mL = RoyiIcholMex(mA, discardThr, 0, shiftVal, maxNumNz); %<! Should match 'IncompleteCholeskyDecompositionMex()'
disp(['Error Norm - ', num2str(norm(mA * vX - (mL * mL.') * vX, 'fro'))]);
mLC = full(mL);
mL = RoyiIcholMex(mA, 0, relThr, shiftVal, maxNumNz); %<! Should match MATLAB's 'ichol()'
disp(['Error Norm - ', num2str(norm(mA * vX - (mL * mL.') * vX, 'fro'))]);
mLD = full(mL);

hF = @() ichol(mA);
hG = @() IncompleteCholeskyDecompositionMex(mA, discardThr, shiftVal, maxNumNz);
hH = @() RoyiIcholMex(mA, discardThr, relThr, shiftVal, maxNumNz);

TimeItMin(hF, 1)
TimeItMin(hG, 1)
TimeItMin(hH, 1)