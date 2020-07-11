% ----------------------------------------------------------------------------------------------- %
% Solve Cholesky Decomposed Linear System Unit Test - 'SolveCholsekyLinearSystem()'
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

funName = 'SolveCholsekyLinearSystem()';


%% Settings

numRows     = 400;
randDensity = 0.75;

discardThr  = 1e-4;
shiftVal    = 0.0;
maxNumNz    = min(5e6, numRows * numRows);


%% Generating Data

mA  = sprandn(numRows, numRows, randDensity);
mA  = (mA.' * mA) + speye(numRows);

vB = randn(numRows, 1);

mL = IncompleteCholeskyDecompositionMex(mA, discardThr, shiftVal, maxNumNz);


%% Analysis

sIchol = struct('type', 'ict', 'droptol', discardThr / 100);

hRunTime = tic();
vX = mL.' \ (mL \ vB);
matlabRunTime = toc(hRunTime);

vXRef = vX;

hRunTime = tic();
vX = SolveCholsekyLinearSystemMex(mL, vB);
dllRunTime = toc(hRunTime);

vE          = vX(:) - vXRef(:);
maxAbsErr   = max(abs(vE));
rmseErr     = sqrt(mean(vE .* vE));

disp([' ']);
disp([funName, ' Unit Test']);
disp(['Max Abs Error    - ', num2str(maxAbsErr)]);
disp(['RMSE Error       - ', num2str(rmseErr)]);
disp(['MATLAB RMSE      - ', num2str(norm(mA * vXRef - vB, 'fro'))]);
disp(['DLL RMSE         - ', num2str(norm(mA * vX - vB, 'fro'))]);
disp(['MATLAB Run Time  - ', num2str(matlabRunTime)]);
disp(['DLL Run Time     - ', num2str(dllRunTime)]);
disp([' ']);


%% Restore Defaults

% set(0, 'DefaultFigureWindowStyle', 'normal');
% set(0, 'DefaultAxesLooseInset', defaultLoosInset);

