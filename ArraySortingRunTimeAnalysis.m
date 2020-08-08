% ----------------------------------------------------------------------------------------------- %
% Incomplete Cholesky Decomposition Run Time Analysis
% Reference:
%   1. fd
% Remarks:
%   1.  Some sources:
%       *   https://github.com/warathul/mergenetsort.
%       *   https://github.com/eggpi/sorting-networks-test.
%       *   https://github.com/983/RadixSort.
%       *   https://github.com/scandum/quadsort.
%       *   https://github.com/scandum/wolfsort.
%       *   https://github.com/swenson/sort/.
%       *   https://www.inf.hs-flensburg.de/lang/algorithmen/sortieren/merge/mergef.htm.
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

DATA_TYPE_UINT32    = 1;
DATA_TYPE_INT32     = 2;
DATA_TYPE_SINGLE    = 3;
DATA_TYPE_DOUBLE    = 4;


%% Settings

vNumElements    = [2:2000]; %<! Goes out of memory for 3000
dataType        = DATA_TYPE_DOUBLE;

cSorter = {@(vA) sort(vA), ...
    @(vA) ArraySortingMex(vA), ...
    @(vA) QuickSortMex(vA), ...
    @(vA) MergeSortMex(vA)};
cSorterString = {['MATLAB Sort'], ['Sorting Network + Manual Quick Sort'], ['LIBC Quick Sort'], ['Merge Sort 4']};

numIterations = 31;


%% Generating Data

numDim      = length(vNumElements);
numSorter   = length(cSorter);
vRunTime    = zeros(numIterations, 1);
mRunTime    = zeros(numDim, numSorter);

for ii = 1:numDim
    disp(['Working on Size: ', num2str(vNumElements(ii),'%05d')]);
    vA = rand(vNumElements(ii), 1) * 2e6;
    vB = zeros(vNumElements(ii), 1);
    vC = zeros(vNumElements(ii), 1);
    vD = zeros(vNumElements(ii), 1);
    
    vB(:) = vA;
    vC(:) = vA;
    vD(:) = vA;
    for jj = 1:numSorter
        disp(['Working on Sorter #', num2str(jj,'%03d'), ' Out of #', num2str(numSorter, '%03d')]);
        if(jj == 1)
            for ll = 1:numIterations
                hRunTime = tic();
                vA = cSorter{jj}(vA);
                vRunTime(ll) = toc(hRunTime);
            end
        elseif(jj == 2)
            for ll = 1:numIterations
                hRunTime = tic();
                cSorter{jj}(vB);
                vRunTime(ll) = toc(hRunTime);
            end
        elseif(jj == 3)
            for ll = 1:numIterations
                hRunTime = tic();
                cSorter{jj}(vC);
                vRunTime(ll) = toc(hRunTime);
            end
        elseif(jj == 4)
            for ll = 1:numIterations
                hRunTime = tic();
                cSorter{jj}(vD);
                vRunTime(ll) = toc(hRunTime);
            end
        end
        disp(['Finished in ', num2str(sum(vRunTime)), ' [Sec]']);
        mRunTime(ii, jj) = min(vRunTime);
    end
end

mRunTime = mRunTime * 1e6;


%% Display Results

figureIdx = figureIdx + 1;

hFigure = figure('Position', [100, 100, 800, 800]);
hAxes = axes();
hLineSeris = plot(vNumElements, mRunTime);
set(hLineSeris, 'LineWidth', 3);
set(get(hAxes, 'Title'), 'String', {['Run Time for Sorting Method']}, 'FontSize', 14);
set(get(hAxes, 'XLabel'), 'String', ['Number of Elements'], 'FontSize', 12);
set(get(hAxes, 'YLabel'), 'String', ['Run Time [Micro Sec]'], 'FontSize', 12);
legend(cSorterString);

if(generateFigures == ON)
    % saveas(hFigure,['Figure', num2str(figureIdx, figureCounterSpec), '.png']);
    print(hFigure, ['Figure', num2str(figureIdx, counterSpec), '.png'], '-dpng', '-r0'); %<! Saves as Screen Resolution
end


%% Auxiliary Functions

function [ mA ] = GenRandPdSparseMat( numRows, randDensity )

% mA  = sprandsym(numRows, randDensity, rand(numRows, 1)) + (5 * speye(numRows)); %<! Very slow!

% Roughly ensuring the diagonal element in each row is lkarger than the
% absolute sum of all other elements in the row
mA  = sprandsym(numRows, randDensity) + (numRows * speye(numRows));


end



%% Restore Defaults

% set(0, 'DefaultFigureWindowStyle', 'normal');
% set(0, 'DefaultAxesLooseInset', defaultLoosInset);

