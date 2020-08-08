% ----------------------------------------------------------------------------------------------- %
% PCGCholesky Unit Test - 'PCGCholesky()'
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

funName = 'PCGCholesky()';


%% Settings

numRows     = 400;
randDensity = 0.75;

discardThr  = 1e-4;
maxNumNz    = min(5e6, numRows * numRows);

absTol  = 1e-16;
relTol  = 1e-18;
maxIter = 100;


%% Generating Data

mA  = sprandn(numRows, numRows, randDensity);
mA  = (mA.' * mA) + speye(numRows);
vB  = randn(numRows, 1);
vX0 = zeros(numRows, 1);

% sIchol = struct('type', 'ict', 'droptol', discardThr / 1000);
% mL = ichol(mA, sIchol);

mL = IncompleteCholeskyDecompositionMex(mA, discardThr, maxNumNz);


%% Analysis

hRunTime = tic();
vX = pcg(mA, vB, 1e-11, maxIter, mL, mL.', vX0);
% vX = pcg(mA, vB, 1e-20, maxIter);
matlabRunTime = toc(hRunTime);

vXRef = vX;

hRunTime = tic();
vX = PCGCholesky(mA, vB, vX0, mL, absTol, relTol, maxIter);
dllRunTime = toc(hRunTime);

vE          = vX(:) - vXRef(:);
maxAbsErr   = max(abs(vE));
rmseErr     = sqrt(mean(vE .* vE));

disp([' ']);
disp([funName, ' Unit Test']);
disp(['Max Abs Error    - ', num2str(maxAbsErr)]);
disp(['RMSE Error       - ', num2str(rmseErr)]);
disp(['MATLAB RMSE      - ', num2str(norm(mA * vXRef - vB))]);
disp(['DLL RMSE         - ', num2str(norm(mA * vX - vB))]);
disp(['MATLAB Run Time  - ', num2str(matlabRunTime)]);
disp(['DLL Run Time     - ', num2str(dllRunTime)]);
disp([' ']);


%% Restore Defaults

% set(0, 'DefaultFigureWindowStyle', 'normal');
% set(0, 'DefaultAxesLooseInset', defaultLoosInset);

