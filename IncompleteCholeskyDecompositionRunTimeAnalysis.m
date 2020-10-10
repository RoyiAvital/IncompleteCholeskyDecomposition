% ----------------------------------------------------------------------------------------------- %
% Incomplete Cholesky Decomposition Run Time Analysis
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

generateFigures = ON;

figureIdx       = 0;
counterSpec     = '%04d';

funName = 'IncompleteCholeskyDecomposition()';

ICHOL_T_ABS_GLOBAL = 1;
ICHOL_T_REL_COLUMN = 2;


%% Settings

vNumRows            = [0100, 0250, 0500, 0750, 1000, 1500, 2000, 2500]; %<! Goes out of memory for 3000
vNumRows            = [0100, 0250, 0375, 0500, 0625, 0750, 0875, 1000]; %<! Goes out of memory for 3000
randDensity         = 0.000005; %<! Multiply it by (8 * vNumRows(end)^4) and make sure it is smaller than memory capacity
randDensityFactor   = 3;

discardThr  = 1e-3;
shiftVal    = 0.001;
maxNumNz    = min((vNumRows(end) * vNumRows(end)) ^ 2, double(intmax('int32')));

% Using numRows ^ 2 as 'gallery('poisson', numRows)' generates matrix of
% size numRows^2 * numRows^2.
hF = @(numRows) GenRandPdSparseMat(numRows ^ 2, randDensity);
hG = @(numRows) GenRandPdSparseMat(numRows ^ 2, randDensityFactor * randDensity); 

cDecomposer = {@(mA) IncompleteCholeskyDecompositionMex(mA, discardThr, shiftVal, maxNumNz, ICHOL_T_ABS_GLOBAL), ...
    @(mA) IncompleteCholeskyDecompositionMex(mA, discardThr, shiftVal, maxNumNz, ICHOL_T_REL_COLUMN), ...
    @(mA) ichol(mA, struct('type', 'nofill', 'droptol', 0, 'michol', 'off', 'shape', 'lower')), ...
    @(mA) ichol(mA, struct('type', 'ict', 'droptol', discardThr, 'michol', 'off', 'shape', 'lower'))};
cDecomposerString = {['MEX ICT Global'], ['MEX ICT Column'], ['MATLAB Zero Fill'], ['MATLAB ICT Column']};
cMatrixGen = {@(numRows) gallery('poisson', numRows), @(numRows) gallery('tridiag', numRows ^ 2), @(numRows) GenWlsMatrix(numRows), hF, hG};
cMatrixType = {['Poisson'], ['Tri Diagonal'], ['Weighted Least Squares (WLS)'], ['Random - ', num2str(randDensity)], ['Random - ', num2str(randDensityFactor * randDensity)]};

numIterations = 5;


%% Generating Data

numDim = length(vNumRows);
numDecomposer = length(cDecomposer);
numMat = length(cMatrixGen);
vRunTime = zeros(numIterations, 1);
mRunTime = zeros(numDim, numDecomposer, numMat);

for ii = 1:numDim
    disp(['Working on Size: ', num2str(vNumRows(ii),'%05d')]);
    for kk = 1:numMat
        mA = cMatrixGen{kk}(vNumRows(ii));
        disp(['Working on Matrix Type: ', cMatrixType{kk}]);
        for jj = 1:numDecomposer
            disp(['Working on Decomposer #', num2str(jj,'%03d'), ' Out of #', num2str(numDecomposer, '%03d')]);
            for ll = 1:numIterations
                hRunTime = tic();
                mL = cDecomposer{jj}(mA);
                vRunTime(ll) = toc(hRunTime);
                mL(1, 1) = 1;
            end
            disp(['Finished in ', num2str(sum(vRunTime)), ' [Sec]']);
            mRunTime(ii, jj, kk) = min(vRunTime);
        end
    end
end


%% Display Results

figureIdx = figureIdx + 1;

hFigure = figure('Position', [100, 100, 800, 1200]);
for ii = 1:numMat
    hAxes = subplot(numMat, 1, ii);
    hLineSeris = plot(vNumRows .^ 2, mRunTime(:, :, ii));
    set(hLineSeris, 'LineWidth', 3);
    set(get(hAxes, 'Title'), 'String', {['Run Time for Matrix Type: ', cMatrixType{ii}]}, 'FontSize', 14);
    set(get(hAxes, 'XLabel'), 'String', ['Number of Rows'], 'FontSize', 12);
    set(get(hAxes, 'YLabel'), 'String', ['Run Time [Sec]'], 'FontSize', 12);
    legend(cDecomposerString);
end

if(generateFigures == ON)
    % saveas(hFigure,['Figure', num2str(figureIdx, figureCounterSpec), '.png']);
    print(hFigure, ['Figure', num2str(figureIdx, counterSpec), '.png'], '-dpng', '-r0'); %<! Saves as Screen Resolution
end


%% Auxiliary Functions

function [ mA ] = GenRandPdSparseMat( numRows, randDensity )

% mA  = sprandsym(numRows, randDensity, rand(numRows, 1)) + (5 * speye(numRows)); %<! Very slow!

% Roughly ensuring the diagonal element in each row is larger than the
% absolute sum of all other elements in the row
mA  = sprandsym(numRows, randDensity) + (numRows * speye(numRows));


end



%% Restore Defaults

% set(0, 'DefaultFigureWindowStyle', 'normal');
% set(0, 'DefaultAxesLooseInset', defaultLoosInset);

