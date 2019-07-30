function lsrPerf = analyzeConcatLsrLog(lg,lidx,filters,cfg,grid)

% lsrPerf = analyzeConcatLsrLog(lg,lidx,filters,cfg)
% lg is concatenated lg over multiple sessions, typically the output of
% concatLogsLsr(). lidx is the master location ID in case there are
% multiple locations in lg (e.g. when analyzing whole grid), cfg is
% analysis configuration structure, filters are filters used to generate
% lgs, for record keeping
% output is structure containing summary stats for multiple indicators,
% such as psychometrics, lgistic regression, motor parameters, view angle

masterTimer = tic;
fprintf('ANALYZING OPTO INACTIVATION...\n')

%% defaults
if nargin < 2; lidx       = [];         end
if nargin < 3; filters    = [];         end
if nargin < 4; cfg        = struct([]); end
if nargin < 5; grid       = [];         end

cfg                       = populateCfg(cfg);

%% get rid of undesired locations and epochs
if ~isempty(lidx); lg     = locationLg(lg,lidx); end
if isfield(filters,'epoch'); epc = filters.epoch; else; epc = []; end
if ~isempty(epc);  lg     = epochLg(lg,epc); end

%% exclude mice with low trial counts if necessary
lg                        = excludeLowTrialCtMice(lg,cfg.minNumTrialsPerMouse);

%% separate control trials with deleted towers if applicable
if sum(strcmpi(lg.laserEpoch,'cueHalf1')) > 0
  whichCueHalf            = 'cueHalf1';
elseif sum(strcmpi(lg.laserEpoch,'cueHalf2')) > 0
  whichCueHalf            = 'cueHalf2';
else
  whichCueHalf            = [];
end
[lg,lghalf,lghalf_vsCtrl] = separateCueHalfCtrls(lg,whichCueHalf);

%% get summary stats (n trials, n mice etc)
lsrPerf.info.grid         = grid;
lsrPerf.filters           = filters;
lsrPerf.cfg               = cfg;
lsrPerf.stats             = getSummaryStats(lg);

%% basic performance indicators (by mouse and overall)
tic; fprintf('\tgetting basic performance indicators...\n')
lsrPerf                   = getBasicPerf(lg,lsrPerf,cfg);
fprintf('\t\tdone after %1.1f min\n',toc/60)

%% view angle
tic; fprintf('\tanalyzing view angles (decoding bootstrap)...')
lsrPerf.viewAngle         = viewAngleDecodeBoot(lg,cfg.posBins,[],cfg.nBootIter);
fprintf(' done after %1.1f min\n',toc/60)

%% motor parameters
tic; fprintf('\tanalyzing unusual motor events...')
lsrPerf                   = getMotorIndicators(lg,lsrPerf,cfg);
fprintf(' done after %1.1f min\n',toc/60)

%% cueHalf controls if applicable
if ~isempty(whichCueHalf)
  thist = tic; fprintf('\tanalyzing cue half control (vs. laser)...\n')
  lsrPerf.cueHalf         = cueHalfCtrl(lghalf,cfg);
  fprintf('\t\tdone after %1.1f min\n',toc(thist)/60)
  
  thist = tic; fprintf('\tanalyzing cue half control (vs. ctrl)...\n')
  lsrPerf.cueHalf_vsCtrl  = cueHalfCtrl(lghalf_vsCtrl,cfg);
  fprintf('\t\tdone after %1.1f min\n',toc(thist)/60)
end

%% stats across mice
tic; fprintf('\tx-mouse stats...')
statsvecs = {'percCorrect','percFinish','bias','speed','excessTravel',...
             'dur','largeViewAngles','earlyTurns','unusualEvents'};
for ii = 1:numel(statsvecs)
  eval(sprintf('[lsrPerf.stats.p_xmouse_%s,lsrPerf.stats.test_xmouse_%s]=calculatePValue(lsrPerf.ctrl.%s_mouse,lsrPerf.lsr.%s_mouse,''classicPaired'',0.05);',...
               statsvecs{ii},statsvecs{ii},statsvecs{ii},statsvecs{ii}));
  if isempty(whichCueHalf); continue; end
  eval(sprintf('[lsrPerf.stats.p_xmouse_vsCueHalf_%s,lsrPerf.stats.test_xmouse_vsCueHalf_%s]=calculatePValue(lsrPerf.cueHalf.ctrl.%s_mouse,lsrPerf.cueHalf.lsr.%s_mouse,''classicPaired'',0.05);',...
               statsvecs{ii},statsvecs{ii},statsvecs{ii},statsvecs{ii}));

end
fprintf(' done after %1.1f sec\n',toc)

%% if this is a no-tower maze, skip logistic regression and psychometrics
if isfield(filters,'mazetype')
  if sum(strcmpi(filters.mazetype,'memGuide_noTowers')) > 0
    fprintf('DONE AFTER %1.1f MIN\n',toc(masterTimer)/60)
    return
  end
end

%% psychometrics
tic; fprintf('\tcomputing pyschometric function...')
lsrPerf.ctrl.psychometric = psychometricFit(lg.choice(~lg.laserON),lg.nCues_RminusL(~lg.laserON),false,cfg.psychBins,true);
lsrPerf.lsr.psychometric  = psychometricFit(lg.choice(lg.laserON),lg.nCues_RminusL(lg.laserON),false,cfg.psychBins,true);
% do by mouse for epoch specific
if ~strcmpi(lg.laserEpoch,'whole') 
  midx = unique(lg.mouseID);
  for iMouse = 1:numel(midx)
    lsrPerf.ctrl.psychometric.mouse(iMouse) = psychometricFit(lg.choice(~lg.laserON & lg.mouseID==midx(iMouse)),  ...
                                                              lg.nCues_RminusL(~lg.laserON & lg.mouseID==midx(iMouse)),false,cfg.psychBins,true);
    lsrPerf.lsr.psychometric.mouse(iMouse)  = psychometricFit(lg.choice(lg.laserON & lg.mouseID==midx(iMouse)),   ...
                                                              lg.nCues_RminusL(lg.laserON & lg.mouseID==midx(iMouse)),false,cfg.psychBins,true);
  end
end
fprintf(' done after %1.1f sec\n',toc)

%% logistic regression
tic; fprintf('\tcomputing and bootstrapping logistic regression...')
lsrPerf                   = bootstrapLogReg(lg,lsrPerf,cfg);
fprintf(' done after %1.1f min\n',toc/60)

%% count total time
fprintf('DONE AFTER %1.1f MIN\n',toc(masterTimer)/60)
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
%% get desired epoch only
function lg = epochLg(lg,epc)

idx         = ~lg.laserON | (lg.laserON & strcmpi(lg.laserEpoch,epc));
vectorNames = fieldnames(lg);
for iVec = 1:numel(vectorNames)
  if isfield(lg,vectorNames{iVec}) && ~strcmpi(vectorNames{iVec},'keyFrameLabels')
    lg.(vectorNames{iVec})(~idx) = [];
  end
end

end

%% ------------------------------------------------------------------------
%% separate tower deletion control trials
function [lg,lghalf,lghalf_vsCtrl] = separateCueHalfCtrls(lg,whichCueHalf)

if isempty(whichCueHalf); lghalf = []; lghalf_vsCtrl = []; return; end

switch whichCueHalf
  case 'cueHalf1'
    isCtrlTrial = cellfun(@(x,y)(sum(x < 100) == 0 & sum(y < 100) == 0),lg.cuePos_R,lg.cuePos_L) & ~lg.laserON;
    
  case 'cueHalf2'
    isCtrlTrial = cellfun(@(x,y)(sum(x > 100) == 0 & sum(y > 100) == 0),lg.cuePos_R,lg.cuePos_L) & ~lg.laserON;

end

vectorNames = fieldnames(lg);
alllg       = lg;
clear lg

for iVec = 1:numel(vectorNames)
  if isfield(alllg,vectorNames{iVec}) && ~strcmpi(vectorNames{iVec},'keyFrameLabels')
    lg.(vectorNames{iVec})            = alllg.(vectorNames{iVec})(~isCtrlTrial);
    lghalf.(vectorNames{iVec})        = alllg.(vectorNames{iVec})(isCtrlTrial | alllg.laserON);
    lghalf_vsCtrl.(vectorNames{iVec}) = alllg.(vectorNames{iVec})(~alllg.laserON);
  end
end

% hack: to use existing code, in lghalf_vsCtrl .ctrl regular trials and
% .laser will be control deleted trials
lghalf_vsCtrl.laserON(isCtrlTrial(~alllg.laserON)) = true;

end

%% ------------------------------------------------------------------------
%% get summary stats
function summaryStats = getSummaryStats(lg)

mouseIDs                        = unique(lg.mouseID);
summaryStats.whichMice          = analysisParams.mice(mouseIDs);
summaryStats.nmice              = numel(mouseIDs);
summaryStats.nsessions          = numel(unique(100.*lg.mouseID + lg.sessionID));
summaryStats.nLsrTrials         = sum(lg.laserON);
summaryStats.nCtrlTrials        = sum(~lg.laserON);
summaryStats.lsrTrialsPerMouse  = zeros(1,summaryStats.nmice);
summaryStats.ctrlTrialsPerMouse = zeros(1,summaryStats.nmice);
for iMouse = 1:summaryStats.nmice
  summaryStats.lsrTrialsPerMouse(iMouse)  = sum(lg.laserON(lg.mouseID == mouseIDs(iMouse)));
  summaryStats.ctrlTrialsPerMouse(iMouse) = sum(lg.mouseID == mouseIDs(iMouse)) - summaryStats.lsrTrialsPerMouse(iMouse);
end
 
end

%% ------------------------------------------------------------------------
%% exclude mice with low trial counts
function lg = excludeLowTrialCtMice(lg,minNT)

delidx = [];
mice   = unique(lg.mouseID);
for mm = 1:numel(mice)
  tidx = find(lg.mouseID(lg.laserON) == mice(mm));
  if numel(tidx) < minNT
    delidx = [delidx find(lg.mouseID == mice(mm))]; 
  end
end

vectorNames = fieldnames(lg);
for iVec = 1:numel(vectorNames)
  if isfield(lg,vectorNames{iVec}) && ~strcmpi(vectorNames{iVec},'keyFrameLabels')
    lg.(vectorNames{iVec})(delidx) = [];
  end
end

end

%% ------------------------------------------------------------------------
%% basic performance indicators
function lsrPerf = getBasicPerf(lg,lsrPerf,cfg)

lsrPerf.stats.stdInterval = normcdf(1, 0, 1) - normcdf(-1, 0, 1);
trialDiffic               = 1 - abs(lg.nCues_RminusL) ./ lg.nCues_total;
isEasy                    = trialDiffic < median(trialDiffic);
isHard                    = trialDiffic > median(trialDiffic);
isLapse                   = abs(lg.nCues_RminusL) >= 10;
% mouse stats
fprintf('\t\tcalculating mouse stats...\n')
mice = unique(lg.mouseID);
for mm = 1:numel(mice)
  % percent completed trials
  lsrPerf.ctrl.percFinish_mouse(mm) = ...
    sum(lg.choice(~lg.laserON & lg.mouseID == mice(mm))~=analysisParams.nilCode)./sum(~lg.laserON & lg.mouseID == mice(mm))*100;
  lsrPerf.lsr.percFinish_mouse(mm)  = ...
    sum(lg.choice(lg.laserON & lg.mouseID == mice(mm))~=analysisParams.nilCode)./sum(lg.laserON & lg.mouseID == mice(mm))*100;

  % percent correct and bias
  [lsrPerf.ctrl.percCorrect_mouse(mm),lsrPerf.ctrl.percCorrect_mouse_binoErr(:,mm)] = ...
    binofit(sum(lg.trialType(~lg.laserON & lg.mouseID == mice(mm)) == lg.choice(~lg.laserON & lg.mouseID == mice(mm))),...
    lsrPerf.stats.ctrlTrialsPerMouse(mm),1-lsrPerf.stats.stdInterval);
  [lsrPerf.lsr.percCorrect_mouse(mm),lsrPerf.lsr.percCorrect_mouse_binoErr(:,mm)] = ...
    binofit(sum(lg.trialType(lg.laserON & lg.mouseID == mice(mm)) == lg.choice(lg.laserON & lg.mouseID == mice(mm))),...
    lsrPerf.stats.lsrTrialsPerMouse(mm),1-lsrPerf.stats.stdInterval);
  
  try
    [lsrPerf.ctrl.percCorrect_easy_mouse(mm),lsrPerf.ctrl.percCorrect_easy_mouse_binoErr(:,mm)] = ...
      binofit(sum(lg.trialType(~lg.laserON & isEasy & lg.mouseID == mice(mm)) == lg.choice(~lg.laserON & isEasy & lg.mouseID == mice(mm))),...
      sum(~lg.laserON & isEasy & lg.mouseID == mice(mm)),1-lsrPerf.stats.stdInterval);
    [lsrPerf.lsr.percCorrect_easy_mouse(mm),lsrPerf.lsr.percCorrect_easy_mouse_binoErr(:,mm)] = ...
      binofit(sum(lg.trialType(lg.laserON & isEasy  & lg.mouseID == mice(mm)) == lg.choice(lg.laserON & isEasy & lg.mouseID == mice(mm))),...
      sum(lg.laserON & isEasy & lg.mouseID == mice(mm)),1-lsrPerf.stats.stdInterval);
  catch
    lsrPerf.ctrl.percCorrect_easy_mouse(mm)           = nan;
    lsrPerf.ctrl.percCorrect_easy_mouse_binoErr(:,mm) = [nan; nan];
    lsrPerf.lsr.percCorrect_easy_mouse(mm)            = nan;
    lsrPerf.lsr.percCorrect_easy_mouse_binoErr(:,mm)  = [nan; nan];
  end
  
  try
    [lsrPerf.ctrl.percCorrect_hard_mouse(mm),lsrPerf.ctrl.percCorrect_hard_mouse_binoErr(:,mm)] = ...
      binofit(sum(lg.trialType(~lg.laserON & isHard & lg.mouseID == mice(mm)) == lg.choice(~lg.laserON & isHard & lg.mouseID == mice(mm))),...
      sum(~lg.laserON & isEasy & lg.mouseID == mice(mm)),1-lsrPerf.stats.stdInterval);
    [lsrPerf.lsr.percCorrect_hard_mouse(mm),lsrPerf.lsr.percCorrect_hard_mouse_binoErr(:,mm)] = ...
      binofit(sum(lg.trialType(lg.laserON & isHard  & lg.mouseID == mice(mm)) == lg.choice(lg.laserON & isHard & lg.mouseID == mice(mm))),...
      sum(lg.laserON & isHard & lg.mouseID == mice(mm)),1-lsrPerf.stats.stdInterval);
  catch
    lsrPerf.ctrl.percCorrect_hard_mouse(mm)           = nan;
    lsrPerf.ctrl.percCorrect_hard_mouse_binoErr(:,mm) = [nan; nan];
    lsrPerf.lsr.percCorrect_hard_mouse(mm)            = nan;
    lsrPerf.lsr.percCorrect_hard_mouse_binoErr(:,mm)  = [nan; nan];
  end
  
  [lsrPerf.ctrl.percCorrect_R_mouse(mm),lsrPerf.ctrl.percCorrect_R_mouse_binoErr(:,mm)] = ...
    binofit(sum(lg.trialType(~lg.laserON & lg.mouseID == mice(mm) & lg.trialType == analysisParams.rightCode) == lg.choice(~lg.laserON & lg.mouseID == mice(mm) & lg.trialType == analysisParams.rightCode)),...
    sum(~lg.laserON & lg.mouseID == mice(mm) & lg.trialType == analysisParams.rightCode),1-lsrPerf.stats.stdInterval);
  [lsrPerf.lsr.percCorrect_R_mouse(mm),lsrPerf.lsr.percCorrect_R_mouse_binoErr(:,mm)] = ...
    binofit(sum(lg.trialType(lg.laserON & lg.mouseID == mice(mm) & lg.trialType == analysisParams.rightCode) == lg.choice(lg.laserON & lg.mouseID == mice(mm) & lg.trialType == analysisParams.rightCode)),...
    sum(lg.laserON & lg.mouseID == mice(mm) & lg.trialType == analysisParams.rightCode),1-lsrPerf.stats.stdInterval);
  
  [lsrPerf.ctrl.percCorrect_L_mouse(mm),lsrPerf.ctrl.percCorrect_L_mouse_binoErr(:,mm)] = ...
    binofit(sum(lg.trialType(~lg.laserON & lg.mouseID == mice(mm) & lg.trialType == analysisParams.leftCode) == lg.choice(~lg.laserON & lg.mouseID == mice(mm) & lg.trialType == analysisParams.leftCode)),...
    sum(~lg.laserON & lg.mouseID == mice(mm) & lg.trialType == analysisParams.leftCode),1-lsrPerf.stats.stdInterval);
  [lsrPerf.lsr.percCorrect_L_mouse(mm),lsrPerf.lsr.percCorrect_L_mouse_binoErr(:,mm)] = ...
    binofit(sum(lg.trialType(lg.laserON & lg.mouseID == mice(mm) & lg.trialType == analysisParams.leftCode) == lg.choice(lg.laserON & lg.mouseID == mice(mm) & lg.trialType == analysisParams.leftCode)),...
    sum(lg.laserON & lg.mouseID == mice(mm) & lg.trialType == analysisParams.leftCode),1-lsrPerf.stats.stdInterval);
  
  [lsrPerf.ctrl.lapse_mouse(mm),lsrPerf.ctrl.lapse_mouse_binoErr(:,mm)] = ...
    binofit(sum(lg.trialType(~lg.laserON & isLapse & lg.mouseID == mice(mm)) ~= lg.choice(~lg.laserON & isLapse & lg.mouseID == mice(mm))),...
    sum(~lg.laserON & isLapse & lg.mouseID == mice(mm)),1-lsrPerf.stats.stdInterval);
  [lsrPerf.lsr.lapse_mouse(mm),lsrPerf.lsr.lapse_mouse_binoErr(:,mm)] = ...
    binofit(sum(lg.trialType(lg.laserON & isLapse & lg.mouseID == mice(mm)) ~= lg.choice(lg.laserON & isLapse & lg.mouseID == mice(mm))),...
    sum(lg.laserON & isLapse & lg.mouseID == mice(mm)),1-lsrPerf.stats.stdInterval);

end

% convert performance to percentage
lsrPerf.ctrl.percCorrect_mouse                = lsrPerf.ctrl.percCorrect_mouse.*100;
lsrPerf.ctrl.percCorrect_mouse_binoErr        = lsrPerf.ctrl.percCorrect_mouse_binoErr.*100;
lsrPerf.lsr.percCorrect_mouse                 = lsrPerf.lsr.percCorrect_mouse.*100;
lsrPerf.lsr.percCorrect_mouse_binoErr         = lsrPerf.lsr.percCorrect_mouse_binoErr.*100;
lsrPerf.lsr.percCorrect_mouse_norm            = (lsrPerf.lsr.percCorrect_mouse - lsrPerf.ctrl.percCorrect_mouse) ./ (lsrPerf.ctrl.percCorrect_mouse - 50);
lsrPerf.lsr.percCorrect_easy_mouse            = lsrPerf.lsr.percCorrect_easy_mouse.*100;
lsrPerf.lsr.percCorrect_easy_mouse_binoErr    = lsrPerf.lsr.percCorrect_easy_mouse_binoErr.*100;
lsrPerf.lsr.percCorrect_hard_mouse            = lsrPerf.lsr.percCorrect_hard_mouse.*100;
lsrPerf.lsr.percCorrect_hard_mouse_binoErr    = lsrPerf.lsr.percCorrect_hard_mouse_binoErr.*100;
lsrPerf.ctrl.percCorrect_R_mouse              = lsrPerf.ctrl.percCorrect_R_mouse.*100;
lsrPerf.ctrl.percCorrect_R_mouse_binoErr      = lsrPerf.ctrl.percCorrect_R_mouse_binoErr.*100;
lsrPerf.lsr.percCorrect_R_mouse               = lsrPerf.lsr.percCorrect_R_mouse.*100;
lsrPerf.lsr.percCorrect_R_mouse_binoErr       = lsrPerf.lsr.percCorrect_L_mouse_binoErr.*100;
lsrPerf.ctrl.percCorrect_L_mouse              = lsrPerf.ctrl.percCorrect_L_mouse.*100;
lsrPerf.ctrl.percCorrect_L_mouse_binoErr      = lsrPerf.ctrl.percCorrect_L_mouse_binoErr.*100;
lsrPerf.lsr.percCorrect_L_mouse               = lsrPerf.lsr.percCorrect_L_mouse.*100;
lsrPerf.lsr.percCorrect_L_mouse_binoErr       = lsrPerf.lsr.percCorrect_L_mouse_binoErr.*100;
lsrPerf.ctrl.lapse_mouse                      = lsrPerf.ctrl.lapse_mouse.*100;
lsrPerf.ctrl.lapse_mouse_binoErr              = lsrPerf.ctrl.lapse_mouse_binoErr.*100;
lsrPerf.lsr.lapse_mouse                       = lsrPerf.lsr.lapse_mouse.*100;
lsrPerf.lsr.lapse_mouse_binoErr               = lsrPerf.lsr.lapse_mouse_binoErr.*100;

% bias (positive value are R)
lsrPerf.ctrl.bias_mouse         = lsrPerf.ctrl.percCorrect_R_mouse - lsrPerf.ctrl.percCorrect_L_mouse;
lsrPerf.ctrl.bias_abs_mouse     = abs(lsrPerf.ctrl.percCorrect_R_mouse - lsrPerf.ctrl.percCorrect_L_mouse);
lsrPerf.lsr.bias_mouse          = lsrPerf.lsr.percCorrect_R_mouse  - lsrPerf.lsr.percCorrect_L_mouse;
lsrPerf.lsr.bias_abs_mouse      = abs(lsrPerf.lsr.percCorrect_R_mouse  - lsrPerf.lsr.percCorrect_L_mouse);

% overall (pooling mice together)
lsrPerf.ctrl.percFinish = sum(lg.choice(~lg.laserON)~=analysisParams.nilCode)./sum(~lg.laserON)*100;
lsrPerf.lsr.percFinish  = sum(lg.choice(lg.laserON)~=analysisParams.nilCode)./sum(lg.laserON)*100;
[~,lsrPerf.ctrl.percFinish_binoErr] = ...
    binofit(sum(lg.choice(~lg.laserON)~=analysisParams.nilCode),lsrPerf.stats.nCtrlTrials,1-lsrPerf.stats.stdInterval);
[~,lsrPerf.lsr.percFinish_binoErr]  = ...
    binofit(sum(lg.choice(lg.laserON)~=analysisParams.nilCode),lsrPerf.stats.nLsrTrials,1-lsrPerf.stats.stdInterval);
  
lsrPerf.ctrl.percCorrect             = ...
  sum(lg.trialType(~lg.laserON) == lg.choice(~lg.laserON))./lsrPerf.stats.nCtrlTrials*100;
lsrPerf.lsr.percCorrect              = ...
  sum(lg.trialType(lg.laserON) == lg.choice(lg.laserON))./lsrPerf.stats.nLsrTrials*100;
lsrPerf.lsr.percCorrect_norm         = ...
  (lsrPerf.lsr.percCorrect - lsrPerf.ctrl.percCorrect) ./ (lsrPerf.ctrl.percCorrect - 50);

[~,lsrPerf.ctrl.percCorrect_binoErr] = ...
    binofit(sum(lg.trialType(~lg.laserON) == lg.choice(~lg.laserON)),lsrPerf.stats.nCtrlTrials,1-lsrPerf.stats.stdInterval);
[~,lsrPerf.lsr.percCorrect_binoErr]  = ...
    binofit(sum(lg.trialType(lg.laserON) == lg.choice(lg.laserON)),lsrPerf.stats.nLsrTrials,1-lsrPerf.stats.stdInterval);

lsrPerf.ctrl.percCorrect_easy        = ...
  sum(lg.trialType(~lg.laserON & isEasy) == lg.choice(~lg.laserON & isEasy))./sum(~lg.laserON & isEasy)*100;
lsrPerf.lsr.percCorrect_easy         = ...
  sum(lg.trialType(lg.laserON & isEasy) == lg.choice(lg.laserON & isEasy))./sum(lg.laserON & isEasy)*100;
[~,lsrPerf.ctrl.percCorrect_easy_binoErr] = ...
    binofit(sum(lg.trialType(~lg.laserON & isEasy) == lg.choice(~lg.laserON & isEasy)),sum(~lg.laserON & isEasy),1-lsrPerf.stats.stdInterval);
[~,lsrPerf.lsr.percCorrect_easy_binoErr]  = ...
    binofit(sum(lg.trialType(lg.laserON & isEasy) == lg.choice(lg.laserON & isEasy)),sum(lg.laserON & isEasy),1-lsrPerf.stats.stdInterval);

lsrPerf.ctrl.percCorrect_hard        = ...
  sum(lg.trialType(~lg.laserON & isHard) == lg.choice(~lg.laserON & isHard))./sum(~lg.laserON & isHard)*100;
lsrPerf.lsr.percCorrect_hard         = ...
  sum(lg.trialType(lg.laserON & isHard) == lg.choice(lg.laserON & isHard))./sum(lg.laserON & isHard)*100;
[~,lsrPerf.ctrl.percCorrect_hard_binoErr] = ...
    binofit(sum(lg.trialType(~lg.laserON & isHard) == lg.choice(~lg.laserON & isHard)),sum(~lg.laserON & isHard),1-lsrPerf.stats.stdInterval);
[~,lsrPerf.lsr.percCorrect_hard_binoErr]  = ...
    binofit(sum(lg.trialType(lg.laserON & isHard) == lg.choice(lg.laserON & isHard)),sum(lg.laserON & isHard),1-lsrPerf.stats.stdInterval);  
  
lsrPerf.ctrl.percCorrect_R             = ...
  sum(lg.trialType(~lg.laserON & lg.trialType == analysisParams.rightCode) == lg.choice(~lg.laserON & lg.trialType == analysisParams.rightCode))./...
  sum(~lg.laserON & lg.trialType == analysisParams.rightCode)*100;
lsrPerf.lsr.percCorrect_R              = ...
  sum(lg.trialType(lg.laserON & lg.trialType == analysisParams.rightCode) == lg.choice(lg.laserON & lg.trialType == analysisParams.rightCode))./...
  sum(lg.laserON & lg.trialType == analysisParams.rightCode)*100;
[~,lsrPerf.ctrl.percCorrect_R_binoErr] = ...
    binofit(sum(lg.trialType(~lg.laserON & lg.trialType == analysisParams.rightCode) == lg.choice(~lg.laserON & lg.trialType == analysisParams.rightCode)),...
    sum(~lg.laserON & lg.trialType == analysisParams.rightCode),1-lsrPerf.stats.stdInterval);
[~,lsrPerf.lsr.percCorrect_R_binoErr]  = ...
    binofit(sum(lg.trialType(lg.laserON & lg.trialType == analysisParams.rightCode) == lg.choice(lg.laserON & lg.trialType == analysisParams.rightCode)),...
    sum(lg.laserON & lg.trialType == analysisParams.rightCode),1-lsrPerf.stats.stdInterval);

lsrPerf.ctrl.percCorrect_L             = ...
  sum(lg.trialType(~lg.laserON & lg.trialType == analysisParams.leftCode) == lg.choice(~lg.laserON & lg.trialType == analysisParams.leftCode))./...
  sum(~lg.laserON & lg.trialType == analysisParams.leftCode)*100;
lsrPerf.lsr.percCorrect_L              = ...
  sum(lg.trialType(lg.laserON & lg.trialType == analysisParams.leftCode) == lg.choice(lg.laserON & lg.trialType == analysisParams.leftCode))./...
  sum(lg.laserON & lg.trialType == analysisParams.leftCode)*100;
[~,lsrPerf.ctrl.percCorrect_L_binoErr] = ...
    binofit(sum(lg.trialType(~lg.laserON & lg.trialType == analysisParams.leftCode) == lg.choice(~lg.laserON & lg.trialType == analysisParams.leftCode)),...
    sum(~lg.laserON & lg.trialType == analysisParams.leftCode),1-lsrPerf.stats.stdInterval);
[~,lsrPerf.lsr.percCorrect_L_binoErr]  = ...
    binofit(sum(lg.trialType(lg.laserON & lg.trialType == analysisParams.leftCode) == lg.choice(lg.laserON & lg.trialType == analysisParams.leftCode)),...
    sum(lg.laserON & lg.trialType == analysisParams.leftCode),1-lsrPerf.stats.stdInterval);

lsrPerf.ctrl.lapse             = sum(lg.trialType(~lg.laserON & isLapse) ~= lg.choice(~lg.laserON & isLapse))./sum(~lg.laserON & isLapse)*100;
lsrPerf.lsr.lapse              = sum(lg.trialType(lg.laserON & isLapse) ~= lg.choice(lg.laserON & isLapse))./sum(lg.laserON & isLapse)*100;
[~,lsrPerf.ctrl.lapse_binoErr] = ...
    binofit(sum(lg.trialType(~lg.laserON & isLapse) ~= lg.choice(~lg.laserON & isLapse)),sum(~lg.laserON & isLapse),1-lsrPerf.stats.stdInterval);
[~,lsrPerf.lsr.lapse_binoErr]  = ...
    binofit(sum(lg.trialType(lg.laserON & isLapse) ~= lg.choice(lg.laserON & isLapse)),sum(lg.laserON & isLapse),1-lsrPerf.stats.stdInterval);

  
lsrPerf.ctrl.bias                      = lsrPerf.ctrl.percCorrect_R - lsrPerf.ctrl.percCorrect_L;
lsrPerf.lsr.bias                       = lsrPerf.lsr.percCorrect_R  - lsrPerf.lsr.percCorrect_L;
lsrPerf.ctrl.bias_abs                  = abs(lsrPerf.ctrl.percCorrect_R - lsrPerf.ctrl.percCorrect_L);
lsrPerf.lsr.bias_abs                   = abs(lsrPerf.lsr.percCorrect_R  - lsrPerf.lsr.percCorrect_L);

% bootstrap
fprintf('\t\tbootstrapping performance...\n')
lsrPerf.stats.bootPerf                 = bootstrapPerf(lg.choice,lg.trialType,lg.laserON,          ...
                                                lg.nCues_RminusL,lg.mouseID,lg.sessionID,          ...
                                                cfg.nBootIter,cfg.drawMiceSessions,cfg.fitInBootstrap,trialDiffic);
end

%% ------------------------------------------------------------------------
%% bootstrap logistic regression for stats 
function lsrPerf = bootstrapLogReg(lg,lsrPerf,cfg)

%% delete aborted trials, excess travel etc
lg                        = cleanupConcatLog(lg,0);
warning('off','all')
%% do overall first
lsrPerf.ctrl.logisticReg  = revCorr_logistic                                    ... 
                            (lg.trialType(~lg.laserON),lg.choice(~lg.laserON),  ...
                             lg.cuePos_R(~lg.laserON),lg.cuePos_L(~lg.laserON), ...
                             [],false,cfg.nBins_logReg,true,cfg.nBootIter_logReg);
lsrPerf.lsr.logisticReg   = revCorr_logistic                                    ... 
                            (lg.trialType(lg.laserON),lg.choice(lg.laserON),    ...
                             lg.cuePos_R(lg.laserON),lg.cuePos_L(lg.laserON),   ...
                             [],false,cfg.nBins_logReg,true,cfg.nBootIter_logReg);

%% bootstrap             
ctrlmean       = lsrPerf.ctrl.logisticReg.values;
lsrmean        = lsrPerf.lsr.logisticReg.values;
diffmean       = lsrmean - ctrlmean;
diffvals       = nan(numel(ctrlmean),cfg.nBootIter);
slopectrl      = nan(1,cfg.nBootIter);
slopelsr       = nan(1,cfg.nBootIter);
dIctrl         = nan(1,cfg.nBootIter);
dIlsr          = nan(1,cfg.nBootIter);
avgWeightCtrl  = mean(ctrlmean);
avgWeightLsr   = mean(lsrmean);
avgctrl        = nan(1,cfg.nBootIter);
avglsr         = nan(1,cfg.nBootIter);

x              = toBinCenters(lsrPerf.ctrl.logisticReg.bins)';
fc             = fit(x,ctrlmean,'poly1');
fl             = fit(x,lsrmean,'poly1');

lsrPerf.ctrl.logisticReg.slope      = fc.p1;
lsrPerf.lsr.logisticReg.slope       = fl.p1;

% rectify at 0.001 for ratio;
ctrlmeanr      = ctrlmean;
lsrmeanr       = lsrmean;
ctrlmeanr(ctrlmeanr<.001) = .001;
lsrmeanr(lsrmeanr<.001)   = .001;
decayIndexCtrl = mean(ctrlmeanr(end-floor(numel(ctrlmeanr)/2)+1:end)+1) / mean(ctrlmeanr(1:floor(numel(ctrlmeanr)/2))+1);
decayIndexLsr  = mean(lsrmeanr(end-floor(numel(ctrlmeanr)/2)+1:end)+1) / mean(lsrmeanr(1:floor(numel(ctrlmeanr)/2))+1);

lsrPerf.ctrl.logisticReg.decayIndex = decayIndexCtrl;
lsrPerf.lsr.logisticReg.decayIndex  = decayIndexLsr;

for iBoot = 1:cfg.nBootIter

  idx         = randsample(numel(lg.choice),numel(lg.choice),true);
  choice      = lg.choice(idx);
  trialType   = lg.trialType(idx);
  cuePos_R    = lg.cuePos_R(idx);
  cuePos_L    = lg.cuePos_L(idx);
  laserON     = lg.laserON(idx);
  
  if sum(laserON) == 0
    thisctrl          = nan(size(ctrlmean));
    thislsr           = nan(size(ctrlmean));
    slopectrl(iBoot)  = nan;
    slopelsr(iBoot)   = nan;
  else
    logReg      = revCorr_logistic(trialType(~laserON),choice(~laserON),  ...
                                   cuePos_R(~laserON),cuePos_L(~laserON), ...
                                   [],false,cfg.nBins_logReg,false);
    thisctrl    = logReg.values;
    logReg      = revCorr_logistic(trialType(laserON),choice(laserON),    ...
                                   cuePos_R(laserON),cuePos_L(laserON),   ...
                                   [],false,cfg.nBins_logReg,false);
    thislsr           = logReg.values;
    fc                = fit(x,thisctrl,'poly1');
    fl                = fit(x,thislsr,'poly1');
    slopectrl(iBoot)  = fc.p1;
    slopelsr(iBoot)   = fl.p1;
  end
  
  diffvals(:,iBoot) = thislsr - thisctrl;
  avgctrl(iBoot)    = mean(thisctrl);
  avglsr(iBoot)     = mean(thislsr);
  
  % rectify at 0.001 for ratio;
  thisctrl(thisctrl<.001) = .001;
  thislsr(thislsr<.001)   = .001;
  dIctrl(iBoot)     = mean(thisctrl(end-floor(numel(ctrlmean)/2)+1:end)+1) / mean(thisctrl(1:floor(numel(ctrlmean)/2))+1);
  dIlsr(iBoot)      = mean(thislsr(end-floor(numel(ctrlmean)/2)+1:end)+1) / mean(thislsr(1:floor(numel(ctrlmean)/2))+1);
  
  

end

lsrPerf.lsr.logisticReg.values_diff = diffmean;
lsrPerf.lsr.logisticReg.sem_diff    = nanstd(diffvals,0,2);
lsrPerf.stats.logisticReg_pvals = zeros(1,numel(ctrlmean));
for iBin = 1:numel(ctrlmean)
  lsrPerf.stats.logisticReg_pvals(iBin) = ...
    sum(sign(diffvals(iBin,:)) ~= sign(diffmean(iBin)))/cfg.nBootIter;
end

lsrPerf.lsr.logisticReg.slope_diff_mean = nanmean(slopelsr-slopectrl);
lsrPerf.lsr.logisticReg.slope_diff_std  = nanstd(slopelsr-slopectrl);
lsrPerf.lsr.logisticReg.slope_p         = nansum(sign(slopelsr-slopectrl) ~= ...
                                          sign(lsrPerf.lsr.logisticReg.slope-lsrPerf.ctrl.logisticReg.slope))/cfg.nBootIter;

lsrPerf.ctrl.logisticReg.decayIndex_std      = nanstd(dIctrl);                                        
lsrPerf.lsr.logisticReg.decayIndex_std       = nanstd(dIlsr);                                        
lsrPerf.lsr.logisticReg.decayIndex_diff_mean = nanmean(dIlsr-dIctrl);
lsrPerf.lsr.logisticReg.decayIndex_diff_std  = nanstd(dIlsr-dIctrl);
lsrPerf.lsr.logisticReg.decayIndex_p         = sum(sign(dIlsr-dIctrl) ~= ...
                                               sign(decayIndexLsr-decayIndexCtrl))/cfg.nBootIter;

lsrPerf.ctrl.logisticReg.avgWeight_mean      = nanmean(avgctrl);  
lsrPerf.ctrl.logisticReg.avgWeight_std       = nanstd(avgctrl);                                             
lsrPerf.lsr.logisticReg.avgWeight_mean       = nanmean(avglsr);  
lsrPerf.lsr.logisticReg.avgWeight_std        = nanstd(avglsr);  
lsrPerf.lsr.logisticReg.avgWeight_diff_mean  = nanmean(avglsr-avgctrl);
lsrPerf.lsr.logisticReg.avgWeight_diff_std   = nanstd(avglsr-avgctrl);
lsrPerf.lsr.logisticReg.avgWeight_p          = nansum(sign(avglsr-avgctrl) ~= ...
                                               sign(avgWeightLsr-avgWeightCtrl))/cfg.nBootIter;
warning('on','all')
end

%% ------------------------------------------------------------------------
%% motor indicators 
function lsrPerf = getMotorIndicators(lg,lsrPerf,cfg)

laserON      = lg.laserON;
speedStem    = lg.speedStem;
excessTravel = lg.excessTravel;
trialType    = lg.trialType;
trialDur     = lg.trialDur;
mouseID      = lg.mouseID;
sessionID    = lg.sessionID;

%% speed 
lsrPerf.ctrl.speed     = nanmean(speedStem(~laserON));
lsrPerf.ctrl.speed_sem = nanstd(speedStem(~laserON))./sqrt(lsrPerf.stats.nCtrlTrials-1);
lsrPerf.lsr.speed      = nanmean(speedStem(laserON));
lsrPerf.lsr.speed_sem  = nanstd(speedStem(laserON))./sqrt(lsrPerf.stats.nLsrTrials-1);
[lsrPerf.stats.p_speed, lsrPerf.stats.test_speed]  ...
                       = calculatePValue(speedStem(laserON),speedStem(~laserON),'classic',0.05);

%% excess travel
lsrPerf.ctrl.excessTravel       = 100.*nanmean(excessTravel(~laserON));
lsrPerf.ctrl.excessTravel_sem   = 100.*nanstd(excessTravel(~laserON))./sqrt(lsrPerf.stats.nCtrlTrials-1);
lsrPerf.ctrl.excessTravel_L     = 100.*nanmean(excessTravel(~laserON & trialType == analysisParams.leftCode));
lsrPerf.ctrl.excessTravel_L_sem = 100.*nanstd(excessTravel(~laserON & trialType == analysisParams.leftCode))./sqrt(sum(~laserON & trialType == analysisParams.leftCode)-1);
lsrPerf.ctrl.excessTravel_R     = 100.*nanmean(excessTravel(~laserON & trialType == analysisParams.rightCode));
lsrPerf.ctrl.excessTravel_R_sem = 100.*nanstd(excessTravel(~laserON & trialType == analysisParams.rightCode))./sqrt(sum(~laserON & trialType == analysisParams.rightCode)-1);
lsrPerf.lsr.excessTravel        = 100.*nanmean(excessTravel(laserON));
lsrPerf.lsr.excessTravel_sem    = 100.*nanstd(excessTravel(laserON))./sqrt(lsrPerf.stats.nLsrTrials-1);
lsrPerf.lsr.excessTravel_L      = 100.*nanmean(excessTravel(laserON & trialType == analysisParams.leftCode));
lsrPerf.lsr.excessTravel_L_sem  = 100.*nanstd(excessTravel(laserON & trialType == analysisParams.leftCode))./sqrt(sum(laserON & trialType == analysisParams.leftCode)-1);
lsrPerf.lsr.excessTravel_R      = 100.*nanmean(excessTravel(laserON & trialType == analysisParams.rightCode));
lsrPerf.lsr.excessTravel_R_sem  = 100.*nanstd(excessTravel(laserON & trialType == analysisParams.rightCode))./sqrt(sum(laserON & trialType == analysisParams.rightCode)-1);
[lsrPerf.stats.p_excessTravel, lsrPerf.stats.test_excessTravel]      ...
                                = calculatePValue(excessTravel(laserON),excessTravel(~laserON),'classic',0.05);
[lsrPerf.stats.p_excessTravel_L, lsrPerf.stats.test_excessTravel_L]  ...
                                = calculatePValue(excessTravel(laserON & trialType == analysisParams.leftCode),excessTravel(~laserON & trialType == analysisParams.leftCode),'classic',0.05);
[lsrPerf.stats.p_excessTravel_R, lsrPerf.stats.test_excessTravel_R]  ...
                                = calculatePValue(excessTravel(laserON & trialType == analysisParams.rightCode),excessTravel(~laserON & trialType == analysisParams.rightCode),'classic',0.05);

%% trial duration
lsrPerf.ctrl.dur       = nanmean(trialDur(~laserON));
lsrPerf.ctrl.dur_sem   = nanstd(trialDur(~laserON))/sqrt(lsrPerf.stats.nCtrlTrials-1);
lsrPerf.lsr.dur        = nanmean(trialDur(laserON));
lsrPerf.lsr.dur_sem    = nanstd(trialDur(laserON))/sqrt(lsrPerf.stats.nLsrTrials-1);
[lsrPerf.stats.p_dur, lsrPerf.stats.test_dur]  ...
                       = calculatePValue(trialDur(laserON),trialDur(~laserON),'classic',0.05);

%% unusual events  
% large view angles
metamouseID                  = -1;
largeAngles                  = unusualMotorEvents(lg,metamouseID,'cueTurnAround');
lsrPerf.ctrl.largeViewAngles = 100 * sum(largeAngles(~laserON)) / sum(~laserON);
lsrPerf.lsr.largeViewAngles  = 100 * sum(largeAngles(laserON)) / sum(laserON);

% early turns
earlyTurns                   = unusualMotorEvents(lg,metamouseID,'earlyTurn');
lsrPerf.ctrl.earlyTurns      = 100 * sum(earlyTurns(~laserON)) / sum(~laserON);
lsrPerf.lsr.earlyTurns       = 100 * sum(earlyTurns(laserON)) / sum(laserON);

% overall unusual events
events                       = unusualMotorEvents(lg,metamouseID,'any');
lsrPerf.ctrl.unusualEvents   = 100 * sum(events(~laserON)) / sum(~laserON);
lsrPerf.lsr.unusualEvents    = 100 * sum(events(laserON)) / sum(laserON);

% bootstrap event indicidence for stats
largeVaSign   = sign(lsrPerf.lsr.largeViewAngles - lsrPerf.ctrl.largeViewAngles);
earlyTurnSign = sign(lsrPerf.lsr.earlyTurns - lsrPerf.ctrl.earlyTurns);
eventsSign    = sign(lsrPerf.lsr.unusualEvents - lsrPerf.ctrl.unusualEvents);

lvactrl       = nan(cfg.nBootIter,1);
etnctrl       = nan(cfg.nBootIter,1);
umectrl       = nan(cfg.nBootIter,1);
lvalsr        = nan(cfg.nBootIter,1);
etnlsr        = nan(cfg.nBootIter,1);
umelsr        = nan(cfg.nBootIter,1);

for iBoot = 1:cfg.nBootIter
  mlist       = unique(lg.mouseID);
  nmice       = length(mlist);
  midx        = mlist(randsample(nmice,nmice,true)); % sample with replacement from mice
  
  % have to loop here to ensure unique sessions for each mouse
  idx = [];
  for mm = 1:nmice
    tidx    = find(mouseID == midx(mm));
    nsess   = length(unique(sessionID(tidx)));
    sidx    = randsample(nsess,nsess,true); % sample with replacement from sessions
    for ss = 1:nsess
      thisl = sum(sessionID(tidx)==sidx(ss));
      idx(end+1:end+thisl) = tidx(sessionID(tidx)==sidx(ss));
    end
  end
  idx          = idx(randsample(numel(idx),numel(idx),true)); % sample with replacement from trials
  iearlyTurns  = earlyTurns(idx);
  ilargeAngles = largeAngles(idx);
  ievents      = events(idx);
  ilaserON     = laserON(idx);
  
  lvactrl(iBoot)   = 100 * (sum(ilargeAngles(~ilaserON))/sum(~ilaserON));
  etnctrl(iBoot)   = 100 * (sum(iearlyTurns(~ilaserON))/sum(~ilaserON));
  umectrl(iBoot)   = 100 * (sum(ievents(~ilaserON))/sum(~ilaserON));
  lvalsr(iBoot)    = 100 * (sum(ilargeAngles(ilaserON))/sum(ilaserON));
  etnlsr(iBoot)    = 100 * (sum(iearlyTurns(ilaserON))/sum(ilaserON));
  umelsr(iBoot)    = 100 * (sum(ievents(ilaserON))/sum(ilaserON));
  
 end

lva   = lvalsr - lvactrl;
etn   = etnlsr - etnctrl;
ume   = umelsr - umectrl;

lsrPerf.stats.bootMotor.p_largeViewAngles = sum(sign(lva) ~= largeVaSign) / cfg.nBootIter;
lsrPerf.stats.bootMotor.p_earlyTurns      = sum(sign(etn) ~= earlyTurnSign) / cfg.nBootIter;
lsrPerf.stats.bootMotor.p_unusualEvents   = sum(sign(ume) ~= eventsSign) / cfg.nBootIter;

lsrPerf.ctrl.largeViewAngles_sd = nanstd(lvactrl);
lsrPerf.ctrl.earlyTurns_sd      = nanstd(etnctrl);
lsrPerf.ctrl.unusualEvents_sd   = nanstd(umectrl);

lsrPerf.lsr.largeViewAngles_sd  = nanstd(lvalsr);
lsrPerf.lsr.earlyTurns_sd       = nanstd(etnlsr);
lsrPerf.lsr.unusualEvents_sd    = nanstd(umelsr);

%% by mouse
mice = unique(mouseID);
for mm = 1:numel(mice)
  
  lsrPerf.ctrl.speed_mouse(mm)           = nanmean(speedStem(~laserON & mouseID==mice(mm)));
  lsrPerf.lsr.speed_mouse(mm)            = nanmean(speedStem(laserON & mouseID==mice(mm)));
  
  lsrPerf.ctrl.excessTravel_mouse(mm)    = 100 * nanmean(excessTravel(~laserON & mouseID==mice(mm)));
  lsrPerf.lsr.excessTravel_mouse(mm)     = 100 * nanmean(excessTravel(laserON & mouseID==mice(mm)));
  
  lsrPerf.ctrl.dur_mouse(mm)             = nanmean(trialDur(~laserON & mouseID==mice(mm)));
  lsrPerf.lsr.dur_mouse(mm)              = nanmean(trialDur(laserON & mouseID==mice(mm)));
  
  largeAngles                            = unusualMotorEvents(lg,mice(mm),'cueTurnAround');
  lsrPerf.ctrl.largeViewAngles_mouse(mm) = 100 * sum(largeAngles(~laserON(mouseID==mice(mm)))) / sum(~laserON(mouseID==mice(mm)));
  lsrPerf.lsr.largeViewAngles_mouse(mm)  = 100 * sum(largeAngles(laserON(mouseID==mice(mm)))) / sum(laserON(mouseID==mice(mm)));

  % early turns
  earlyTurns                             = unusualMotorEvents(lg,mice(mm),'earlyTurn');
  lsrPerf.ctrl.earlyTurns_mouse(mm)      = 100 * sum(earlyTurns(~laserON(mouseID==mice(mm)))) / sum(~laserON(mouseID==mice(mm)));
  lsrPerf.lsr.earlyTurns_mouse(mm)       = 100 * sum(earlyTurns(laserON(mouseID==mice(mm)))) / sum(laserON(mouseID==mice(mm)));

  % overall unusual events
  events                                 = unusualMotorEvents(lg,mice(mm),'any');
  lsrPerf.ctrl.unusualEvents_mouse(mm)   = 100 * sum(events(~laserON(mouseID==mice(mm)))) / sum(~laserON(mouseID==mice(mm)));
  lsrPerf.lsr.unusualEvents_mouse(mm)    = 100 * sum(events(laserON(mouseID==mice(mm)))) / sum(laserON(mouseID==mice(mm)));
end

end

%% ------------------------------------------------------------------------
%% cue half control
function cueHalf = cueHalfCtrl(lg,cfg)

%% stats
cueHalf.stats             = getSummaryStats(lg);

%% psychometric
cueHalf.psychometric      = ...
  psychometricFit(lg.choice(~lg.laserON),lg.nCues_RminusL(~lg.laserON),false,cfg.psychBins_cueHalf,true);

%% basic performance indicators (by mouse and overall)
tic; fprintf('\t\tgetting basic performance indicators...\n')
cueHalf                   = getBasicPerf(lg,cueHalf,cfg);
fprintf('\t\tdone after %1.1f min\n',toc/60)

%% view angle
tic; fprintf('\t\tanalyzing view angles (decoding bootstrap)...')
cueHalf.viewAngle         = viewAngleDecodeBoot(lg,cfg.posBins,[],cfg.nBootIter);
fprintf(' done after %1.1f min\n',toc/60)

%% motor parameters
tic; fprintf('\t\tanalyzing unusual motor events...')
cueHalf                   = getMotorIndicators(lg,cueHalf,cfg);
fprintf(' done after %1.1f min\n',toc/60)

end

%% ------------------------------------------------------------------------
%% default config
function cfg = populateCfg(cfg)

if ~isfield(cfg,'minNumTrialsPerMouse')
  cfg(1).minNumTrialsPerMouse = 1;
end
if ~isfield(cfg,'nBootIter')
  cfg(1).nBootIter            = 10000;
end
if ~isfield(cfg,'drawMiceSessions')
  cfg(1).drawMiceSessions     = true;
end
if ~isfield(cfg,'fitInBootstrap')
  cfg(1).fitInBootstrap       = true;
end
if ~isfield(cfg,'psychBins')
  cfg(1).psychBins            = -15:5:15;
end
if ~isfield(cfg,'psychBins_cueHalf')
  cfg(1).psychBins_cueHalf    = -9:3:9;
end
if ~isfield(cfg,'nBins_logReg')
  cfg(1).nBins_logReg         = 2;
end
if ~isfield(cfg,'nBootIter_logReg')
  cfg(1).nBootIter_logReg     = 100;
end
if ~isfield(cfg,'posBins')
  cfg(1).posBins              = 0:5:300;
end

end
