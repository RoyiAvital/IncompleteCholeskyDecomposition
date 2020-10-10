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

ICHOL_T_ABS_GLOBAL = 1;
ICHOL_T_REL_COLUMN = 2; %<! Should match MATLAB


%% Settings

numRows     = 500;
randDensity = 0.75;

discardThr  = 1e-4;
shiftVal    = 0.0;
maxNumNz    = min(5e6, numRows ^ 4);
icholType   = ICHOL_T_REL_COLUMN;


%% Generating Data

mA  = sprandn(numRows, numRows, randDensity);
mA  = (mA.' * mA) + speye(numRows);
% mA  = gallery('poisson', numRows);
% mA  = gallery('tridiag', numRows ^ 2);
% mA  = GenWlsMatrix(numRows);


%% Analysis

sIchol = struct('type', 'ict', 'droptol', discardThr, 'michol', 'off');

hRunTime = tic();
mL = ichol(mA, sIchol);
matlabRunTime = toc(hRunTime);

mLRef = mL;

hRunTime = tic();
mL = IncompleteCholeskyDecompositionMex(mA, discardThr, shiftVal, maxNumNz, icholType);
dllRunTime = toc(hRunTime);

vE          = mL(:) - mLRef(:);
maxAbsErr   = max(abs(vE));
rmseErr     = sqrt(mean(vE .* vE));

disp([' ']);
disp([funName, ' Unit Test']);
disp(['Max Abs Error    - ', num2str(maxAbsErr)]);
disp(['RMSE Error       - ', num2str(rmseErr)]);
disp(['MATLAB RMSE      - ', num2str(norm(mA - (mLRef * mLRef.'), 'fro'))]);
disp(['MEX RMSE         - ', num2str(norm(mA - (mL * mL.'), 'fro'))]);
disp(['MATLAB Run Time  - ', num2str(matlabRunTime)]);
disp(['MEX Run Time     - ', num2str(dllRunTime)]);
disp([' ']);


%% Restore Defaults

% set(0, 'DefaultFigureWindowStyle', 'normal');
% set(0, 'DefaultAxesLooseInset', defaultLoosInset);

