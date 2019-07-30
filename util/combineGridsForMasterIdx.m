function lidx = combineGridsForMasterIdx(gridList)

% lidx = combineGridsForMasterIdx(gridList)
% return list of master location indices for single grid or combinations
% thereof. 'any' will return all single-location grids

if nargin < 1
  gridList = 'any';
end

if ~iscell(gridList) && strcmpi(gridList,'any')
  gridList = {'anteriorGrid.mat','reducedGridPosterior.mat','posteriorGrid.mat',...
    'anteriorExtraGrid.mat','anteriorLessExtraGrid.mat','anteriorMostGrid.mat'};
elseif ~iscell(gridList) && strcmpi(gridList,'anyBilateral')
  gridList = {'fullGridBilateral.mat','fullLessGridBilateral.mat'};
elseif ~iscell(gridList) && strcmpi(gridList,'veryReducedGrid')
  gridList = {'veryReducedGrid.mat';'RSCunilateralGrid.mat'};
elseif ~iscell(gridList) && strcmpi(gridList,'veryReducedGridBilateral')
  gridList = {'veryReducedGridBilateral.mat';'ML10_AP30_bilateralGrid.mat';...
    'ML5_AP-25_bilateralGrid.mat'; 'PPCHarveyBilateralGrid.mat';...
    'veryReducedGridBilateral.mat'};
elseif ~iscell(gridList) && strcmpi(gridList,'veryReducedGridBilateral_v2')
  gridList = {'veryReducedGridBilateral.mat';'ML10_AP30_bilateralGrid.mat';...
    'ML5_AP-25_bilateralGrid.mat'; 'PPCHarveyBilateralGrid.mat';...
    'veryReducedGridBilateral.mat'};
end

if iscell(gridList) && strcmpi(gridList{1},'anteriorGrid.mat')
  lidx = [];
  for ii = 1:length(gridList)
    load([analysisParams.gridpath gridList{ii}],'grid')
    lidx = [lidx; getMasterLocationIdx(grid)];
  end
  lidx = unique(lidx);
elseif iscell(gridList) && ...
    (strcmpi(gridList{1},'fullGridBilateral.mat')       || ...
    strcmpi(gridList{1},'veryReducedGrid.mat')          || ...
    strcmpi(gridList{1},'veryReducedGridBilateral.mat') || ...
    strcmpi(gridList{1},'veryReducedGridBilateral_v2.mat'))
  load([analysisParams.gridpath gridList{1}],'grid')
  lidx = getMasterLocationIdx(grid);
else
  load([analysisParams.gridpath gridList],'grid')
  lidx = getMasterLocationIdx(grid);
end