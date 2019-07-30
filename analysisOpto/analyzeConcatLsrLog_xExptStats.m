function analyzeConcatLsrLog_xExptStats(spockFlag,updateDataset,niter)

% analyzeConcatLsrLog_xExptStats(spockFlag,updateDataset,niter)
% compares significance across experiments

%% initialize
if nargin < 1; spockFlag     = true;   end
if nargin < 2; updateDataset = false;  end
if nargin < 3; niter         = 10000;  end


poolobj = gcp('nocreate');
if isempty(poolobj); 
  if spockFlag
    poolobj = parpool;
  else
    poolobj = parpool(4);
  end
end
ap      = analysisParams;

try
%% update logs, save data diagnosis
if updateDataset; summarizeLsrExpt_batch; end

warning('off','all'); 

%% whole trial, wt vs 2mW vs 6 mW
thist = tic;
fprintf('-------------------\n\n')

% generate logs
filters.mazetype                = {'last';'last-1'};
filters.perfTh                  = 0.6;
filters.deleteNonMainMazeTrials = true;
filters.epoch                   = 'whole';
filters.genotype                = 'vgat';
filters.gridLabel               = {'fullGridBilateral.mat';      ...
                                   'fullLessGridBilateral.mat';  ...     
                                   'posteriorGridBilateral.mat'; ...
                                   'anteriorExtraGridBilateral.mat'};

lg                  = concatLogsLsr([],filters,spockFlag);

filters.powerLogic  = '>=';
filters.power       = 6;
lg6mW               = concatLogsLsr([],filters,spockFlag);

filters.mazetype    = 'visualGuide';
filters.perfTh      = 0.8;
lgt4                = concatLogsLsr([],filters,spockFlag);

filters.mazetype    = {'last';'last-1'};
filters.powerLogic  = '==';
filters.power       = 2;
filters.perfTh      = 0.6;
lg2mW               = concatLogsLsr([],filters,spockFlag);

filters             = rmfield(filters,'powerLogic');
filters             = rmfield(filters,'power');
filters.genotype    = 'wt';
lgwt                = concatLogsLsr([],filters,spockFlag);


data                = load('fullGridBilateral.mat','grid');
grid                = data.grid;
locIdx              = getMasterLocationIdx(grid);

% loop over locations
parfor iLoc = 1:numel(locIdx)
  sublg                 = locationLg(lg,locIdx(iLoc));
  sublg6mW              = locationLg(lg6mW,locIdx(iLoc));
  sublg2mW              = locationLg(lg2mW,locIdx(iLoc));
  sublgwt               = locationLg(lgwt,locIdx(iLoc));
  sublgt4               = locationLg(lgt4,locIdx(iLoc));
  pvals_VGvsCtrl(iLoc,:)= compareEffectSize(sublgwt, sublg, niter, 'onesided');
  pvals_2vs6(iLoc,:)    = compareEffectSize(sublg2mW, sublg6mW, niter, 'onesided');
  pvals_2vsctrl(iLoc,:) = compareEffectSize(sublg2mW, sublgwt,  niter, 'onesided');
  pvals_6vsctrl(iLoc,:) = compareEffectSize(sublg6mW, sublgwt,  niter, 'onesided');
  pvals_t4vst11(iLoc,:) = compareEffectSize(sublg6mW, sublgt4,  niter, 'onesided');
end

% reorganize data structure and FDR correct
pvals.alpha         = 0.05;
pvals.test          = 'onesided';
pvals.grid          = grid;
pvals.vals_2vs6     = pvals_2vs6;
[pvals.isSig_2vs6, pvals.alpha_correct_2vs6] ...
                    = FDR(pvals_2vs6, pvals.alpha);
pvals.vals_2vsCtrl  = pvals_2vsctrl;
[pvals.isSig_2vsCtrl, pvals.alpha_correct_2vsCtrl] ...
                    = FDR(pvals_2vsctrl, pvals.alpha);
pvals.vals_6vsCtrl  = pvals_6vsctrl;
[pvals.isSig_6vsCtrl, pvals.alpha_correct_6vsCtrl] ...
                    = FDR(pvals_6vsctrl, pvals.alpha);
pvals.vals_accumVsVisGuide  = pvals_t4vst11;
[pvals.isSig_accumVsVisGuide, pvals.alpha_accumVsVisGuide] ...
                            = FDR(pvals_t4vst11, pvals.alpha);
pvals.vals_vgatvsCtrl  = pvals_VGvsCtrl;
[pvals.isSig_vgatvsCtrl, pvals.alpha_vgatvsCtrl] ...
                            = FDR(pvals_VGvsCtrl, pvals.alpha);
% save
rootdir         = ap.getRootDir(spockFlag);
fn              = 'fullGrid_wholeTrial_xExptComp';

save([rootdir fn],'pvals')

fprintf('\nANALYZED %s in %1.1f min\n',upper(fn),toc(thist)/60)
clear filters pvals

%% reduced grid, x-epoch and x-loc
thist = tic;
fprintf('-------------------\n\n')

% generate logs
filters.mazetype                = {'last';'last-1'};
filters.perfTh                  = 0.6;
filters.deleteNonMainMazeTrials = true;
filters.powerLogic              = '>=';
filters.power                   = 6;
filters.genotype                = 'vgat';
filters.gridLabel               = {'veryReducedGridBilateral.mat';    ... 
                                   'veryReducedGridBilateral_v2.mat'; ...
                                   'veryReducedGridBilateral_v3.mat'; ...
                                   'ML10_AP25_bilateralGrid.mat';     ...
                                   'ML10_AP30_bilateralGrid.mat';     ...
                                   'PPCHarveyBilateralGrid.mat';      ...
                                   'ML5_AP-25_bilateralGrid.mat';     ...
                                   };

% allow for slightly different locations in case of V1, aM2
loclbls   = {'V1','PPC','RSC','mM2','aM2'};
coords{1} = {[-2.5 -3.5; 2.5 -3.5];[-3 -3.5; 3 -3.5]};
coords{2} = {[-1.75 -2; 1.75 -2]};
coords{3} = {[-.5 -2.5; .5 -2.5]};
coords{4} = {[-.5 0; .5 0]};
coords{5} = {[-1 3; 1 3];[-1 2.5; 1 2.5]};
data      = load('veryReducedGridBilateral_v2.mat','grid');
grid      = data.grid([4 2 5 3 1]);

filters.epoch       = 'cueHalf1';
lgch1               = concatLogsLsr([],filters,spockFlag);

filters.epoch       = 'cueHalf2';
lgch2               = concatLogsLsr([],filters,spockFlag);

filters.epoch       = 'mem';
lgmem               = concatLogsLsr([],filters,spockFlag);

% loop over locations for x-epoch comps
parfor iLoc = 1:numel(coords)
  locIdx                 = getMasterLocationIdx(coords{iLoc});
  sublgch1               = locationLg(lgch1,locIdx);
  sublgch2               = locationLg(lgch2,locIdx);
  sublgmem               = locationLg(lgmem,locIdx);
  pvals_ch1VSch2(iLoc,:) = compareEffectSize(sublgch1, sublgch2, niter);
  pvals_ch1VSdel(iLoc,:) = compareEffectSize(sublgch1, sublgmem, niter);
  pvals_ch2VSdel(iLoc,:) = compareEffectSize(sublgch2, sublgmem, niter);
end

% reorganize data structure and FDR correct
pvals.alpha          =  analysisParams.bootalpha;
pvals.test           = 'onesided';
pvals.grid           = grid;
pvals.loclbls        = loclbls;
pvals.vals_ch1VSch2  = pvals_ch1VSch2;
[pvals.isSig_ch1VSch2, pvals.alpha_correct_ch1VSch2] ...
                     = FDR(pvals_ch1VSch2, analysisParams.bootalpha);
pvals.vals_ch1VSdel  = pvals_ch1VSdel;
[pvals.isSig_ch1VSdel, pvals.alpha_correct_ch1VSdel] ...
                     = FDR(pvals_ch1VSdel, analysisParams.bootalpha);
pvals.vals_ch2VSdel  = pvals_ch2VSdel;
[pvals.isSig_ch2VSdel, pvals.alpha_correct_ch2VSdel] ...
                     = FDR(pvals_ch2VSdel, analysisParams.bootalpha);

% loop over locations for x-location, within epoch comps
pvals_ch1_xLoc = [];
pvals_ch2_xLoc = [];
pvals_mem_xLoc = [];
pairs          = [];
for iLoc = 1:numel(coords)
  locIdx                 = getMasterLocationIdx(coords{iLoc});
  sublgch1_1             = locationLg(lgch1,locIdx);
  sublgch2_1             = locationLg(lgch2,locIdx);
  sublgmem_1             = locationLg(lgmem,locIdx);
  for iLoc2 = 1:numel(coords)
    if iLoc2 <= iLoc; continue; end
    locIdx                  = getMasterLocationIdx(coords{iLoc});
    sublgch1_2              = locationLg(lgch1,locIdx);
    sublgch2_2              = locationLg(lgch2,locIdx);
    sublgmem_2              = locationLg(lgmem,locIdx);
    pvals_ch1_xLoc(end+1,:) = compareEffectSize(sublgch1_1, sublgch1_2, niter);
    pvals_ch2_xLoc(end+1,:) = compareEffectSize(sublgch2_1, sublgch2_2, niter);
    pvals_mem_xLoc(end+1,:) = compareEffectSize(sublgmem_1, sublgmem_2, niter);
    pairs(end+1,:)          = [iLoc iLoc2];
  end
end

% reorganize data structure and FDR correct
pvals.vals_ch1_xLoc  = pvals_ch1_xLoc;
vec                  = pvals_ch1_xLoc(:);
[pvals.isSig_ch1_xLoc, pvals.alpha_correct_ch1_xLoc] ...
                     = FDR(vec(~isnan(vec)), analysisParams.bootalpha);
pvals.vals_ch2_xLoc  = pvals_ch2_xLoc;
vec                  = pvals_ch2_xLoc(:);
[pvals.isSig_ch2_xLoc, pvals.alpha_correct_ch2_xLoc] ...
                     = FDR(vec(~isnan(vec)), analysisParams.bootalpha);
pvals.vals_del_xLoc  = pvals_mem_xLoc;
vec                  = pvals_mem_xLoc(:);
[pvals.isSig_ch2VSdel, pvals.alpha_correct_ch2VSdel] ...
                     = FDR(vec(~isnan(vec)), analysisParams.bootalpha);
                   
% save
rootdir         = ap.getRootDir(spockFlag);
fn              = 'reducedGrid_subTrial_xExptComp';

save([rootdir fn],'pvals')

fprintf('\nANALYZED %s in %1.1f min\n',upper(fn),toc(thist)/60)

catch ME
  displayException(ME)
end
end

%% ------------------------------------------------------------------------
%% get desired location only
function lg = locationLg(lg,lidx)

if iscell(lidx)
  posidx = [];
  for iLoc = 1:numel(lidx)
    posidx(end+1,:) = cellfun(@(x)(isempty(x) || sum(x==lidx{iLoc})==numel(lidx{iLoc})),lg.galvoPosIdxMaster);
  end
  if iLoc > 1; posidx = sum(posidx) > 0; else posidx = posidx > 0; end
else
  posidx = lg.galvoPosIdxMaster == 0;
  for iLoc = 1:numel(lidx)
    posidx(end+1,:) = lg.galvoPosIdxMaster == lidx(iLoc);
  end
  posidx = sum(posidx) > 0;
end

vectorNames = fieldnames(lg);
for iVec = 1:numel(vectorNames)
  if isfield(lg,vectorNames{iVec}) && ~strcmpi(vectorNames{iVec},'keyFrameLabels')
    lg.(vectorNames{iVec})(~posidx) = [];
  end
end
end

%% ------------------------------------------------------------------------
%% get desired location only
function pval = compareEffectSize(lg1,lg2,niter,type)

if nargin < 3; niter = 10000;      end
if nargin < 4; type  = 'twosided'; end % onesided tests lg2 < lg1

% avg effect size
es1 = (sum(lg1.trialType(lg1.laserON) == lg1.choice(lg1.laserON)) / sum(lg1.laserON)) - ...
      (sum(lg1.trialType(~lg1.laserON) == lg1.choice(~lg1.laserON)) / sum(~lg1.laserON));
es2 = (sum(lg2.trialType(lg2.laserON) == lg2.choice(lg2.laserON)) / sum(lg2.laserON)) - ...
      (sum(lg2.trialType(~lg2.laserON) == lg2.choice(~lg2.laserON)) / sum(~lg2.laserON));
    
% generic bootstrapping (two-tailed)
bootvec  = zeros(1,niter);
ttAll    = [lg1.trialType'; lg2.trialType'];
chAll    = [lg1.choice'; lg2.choice'];
lsAll    = [lg1.laserON'; lg2.laserON'];
vecTags  = [true(numel(lg1.trialType),1); false(numel(lg2.trialType),1)];

for iBoot = 1:niter
  idx            = randsample(length(ttAll),length(ttAll),true); % sample with replacement
  ittAll         = ttAll(idx);
  ichAll         = chAll(idx);
  ilsAll         = lsAll(idx);
  ivecTags       = vecTags(idx);
  ie1            = (sum(ittAll(ilsAll & ivecTags) == ichAll((ilsAll & ivecTags))) / sum(ilsAll & ivecTags)) - ...
                   (sum(ittAll(~ilsAll & ivecTags) == ichAll((~ilsAll & ivecTags))) / sum(~ilsAll & ivecTags));
  ie2            = (sum(ittAll(ilsAll & ~ivecTags) == ichAll((ilsAll & ~ivecTags))) / sum(ilsAll & ~ivecTags)) - ...
                   (sum(ittAll(~ilsAll & ~ivecTags) == ichAll((~ilsAll & ~ivecTags))) / sum(~ilsAll & ~ivecTags));
  bootvec(iBoot) = ie2 - ie1;
end

switch type
  case 'twosided'
    pval = sum(sign(bootvec) ~= sign(es2 - es1)) / niter;
  case 'onesided'
    pval = sum(sign(bootvec) > 0) / niter;
end
end