% ----------------------------------------------------------------------------------------------- %
% Array Sorting Unit Test - 'ArraySorting()'
% Reference:
%   1. fd
% Remarks:
%   1.  Working on Float (Single).
% TODO:
%   1.  A
%   Release Notes:
%   -   1.0.000     13/07/2020  Royi Avital
%       *   First release version.
% ----------------------------------------------------------------------------------------------- %

%% Setting Environment Parameters

subStreamNumberDefault = 192;
run('InitScript.m');

funName = 'ArraySorting()';


%% Settings

numElemnts = 128;


%% Generating Data

vA = randn(numElemnts, 1);
vB = zeros(numElemnts, 1);
vC = zeros(numElemnts, 1);
vD = zeros(numElemnts, 1);
vZ = zeros(numElemnts, 1);

vB(:) = vA; %<! To prevent copy by Reference
vC(:) = vA; %<! To prevent copy by Reference
vD(:) = vA; %<! To prevent copy by Reference
vZ(:) = vA; %<! To prevent copy by Reference


%% Analysis

hRunTime = tic();
vB(:) = sort(vB);
matlabRunTime = toc(hRunTime);

hRunTime = tic();
ArraySortingMex(vC);
netRunTime = toc(hRunTime);

vE          = vB(:) - vC(:);
maxAbsErr   = max(abs(vE));
rmseErr     = sqrt(mean(vE .* vE));

disp([' ']);
disp([funName, ' Unit Test']);
disp(['Max Abs Error    - ', num2str(maxAbsErr)]);
disp(['RMSE Error       - ', num2str(rmseErr)]);
disp(['MATLAB Run Time  - ', num2str(matlabRunTime)]);
disp(['Net Run Time     - ', num2str(netRunTime)]);
disp([' ']);

hRunTime    = tic();
QuickSortMex(vD);
quickRunTime = toc(hRunTime);

vE          = vB(:) - vD(:);
maxAbsErr   = max(abs(vE));
rmseErr     = sqrt(mean(vE .* vE));

disp([' ']);
disp([funName, ' Unit Test']);
disp(['Max Abs Error    - ', num2str(maxAbsErr)]);
disp(['RMSE Error       - ', num2str(rmseErr)]);
disp(['MATLAB Run Time  - ', num2str(matlabRunTime)]);
disp(['Quick Run Time   - ', num2str(quickRunTime)]);
disp([' ']);

hRunTime    = tic();
MergeSortMex(vZ);
quickRunTime = toc(hRunTime);

vE          = vB(:) - vZ(:);
maxAbsErr   = max(abs(vE));
rmseErr     = sqrt(mean(vE .* vE));

disp([' ']);
disp([funName, ' Unit Test']);
disp(['Max Abs Error    - ', num2str(maxAbsErr)]);
disp(['RMSE Error       - ', num2str(rmseErr)]);
disp(['MATLAB Run Time  - ', num2str(matlabRunTime)]);
disp(['Merge Run Time   - ', num2str(quickRunTime)]);
disp([' ']);


%% Restore Defaults

% set(0, 'DefaultFigureWindowStyle', 'normal');
% set(0, 'DefaultAxesLooseInset', defaultLoosInset);

