% ----------------------------------------------------------------------------------------------- %
% Incomplete Cholesky Decomposition Unit Test - 'IncompleteCholeskyDecomposition()'
% Reference:
%   1. fd
% Remarks:
%   1.  Working on Float (Single).
% TODO:
%   1.  A
%   Release Notes:
%   -   1.0.000     08/07/2020  Royi Avital
%       *   First release version.
% ----------------------------------------------------------------------------------------------- %

%% Setting Environment Parameters

subStreamNumberDefault = 192;
run('InitScript.m');

funName = 'IncompleteCholeskyDecomposition()';


%% Settings

numRows     = 400;
randDensity = 0.75;

discardThr  = 1e-4;
shiftVal    = 0.0;
maxNumNz    = min(5e6, numRows * numRows);


%% Generating Data

mA  = sprandn(numRows, numRows, randDensity);
mA  = (mA.' * mA) + speye(numRows);


%% Analysis

sIchol = struct('type', 'ict', 'droptol', discardThr / 100);

hRunTime = tic();
mL = ichol(mA, sIchol);
matlabRunTime = toc(hRunTime);

mLRef = mL;

hRunTime = tic();
mL = IncompleteCholeskyDecompositionMex(mA, discardThr, shiftVal, maxNumNz);
dllRunTime = toc(hRunTime);

vE          = mL(:) - mLRef(:);
maxAbsErr   = max(abs(vE));
rmseErr     = sqrt(mean(vE .* vE));

disp([' ']);
disp([funName, ' Unit Test']);
disp(['Max Abs Error    - ', num2str(maxAbsErr)]);
disp(['RMSE Error       - ', num2str(rmseErr)]);
disp(['MATLAB RMSE      - ', num2str(norm(mA - (mLRef * mLRef.'), 'fro'))]);
disp(['DLL RMSE         - ', num2str(norm(mA - (mL * mL.'), 'fro'))]);
disp(['MATLAB Run Time  - ', num2str(matlabRunTime)]);
disp(['DLL Run Time     - ', num2str(dllRunTime)]);
disp([' ']);


%% Restore Defaults

% set(0, 'DefaultFigureWindowStyle', 'normal');
% set(0, 'DefaultAxesLooseInset', defaultLoosInset);

