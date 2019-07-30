function varPower = analyzeVarPower(spockFlag)

% varPower = analyzeVarPower(spockFlag)
% analyzes experiments varying power over V1 during visually guided task

if nargin < 1; spockFlag = true; end

%% record some variables
varPower.vgat.mice = {'vg3';'vg4'};
varPower.ctrl.mice = {'wt6';'wt11';'wt12'};
varPower.powers    = [.25 .5 1 2 4 8 12];

varPower.stats.vgat.nmice = numel(varPower.vgat.mice);
varPower.stats.ctrl.nmice = numel(varPower.ctrl.mice);

%% concat logs
filters.gridLabel   = 'V1unilateralSingleGrid.mat';
filters.genotype    = 'vgat';
filters.mazes       = 4;
filters.mazeLogic   = '=='; % '==', '<=', '>=' string w/logical operator
lg                  = concatLogsLsr([],filters,spockFlag);

filters.genotype    = 'wt';
lg_wt               = concatLogsLsr([],filters,spockFlag);

load(filters.gridLabel,'grid')
lidx = getMasterLocationIdx(grid);

%% analyze just % correct and excess travel, lump into ipsi/contra, vgat
stdInterval = normcdf(1, 0, 1) - normcdf(-1, 0, 1);
for iP = 1:numel(varPower.powers)
  
  %% vgat
  choice       = lg.choice(lg.power == varPower.powers(iP));
  trialType    = lg.trialType(lg.power == varPower.powers(iP));
  laserON      = lg.laserON(lg.power == varPower.powers(iP));
  excessTravel = lg.excessTravel(lg.power == varPower.powers(iP));
  locIdx       = cell2mat(lg.galvoPosIdxMaster(lg.power == varPower.powers(iP)));
  
  varPower.stats.vgat.nTrials_ctrl(iP) = sum(~lg.laserON);
  varPower.stats.vgat.nTrials_lsr(iP)  = sum(laserON);
  
  varPower.vgat.percCorrect_ctrl(iP)   = 100 * sum(lg.trialType(~lg.laserON) == lg.choice(~lg.laserON))./sum(~lg.laserON);
  ipsiIdx                              = laserON & ((trialType == analysisParams.leftCode & locIdx == lidx(1)) | ...
                                                    (trialType == analysisParams.rightCode & locIdx == lidx(2)));
  varPower.vgat.percCorrect_ipsi(iP)   = 100 * sum(trialType(ipsiIdx) == choice(ipsiIdx)) / sum(ipsiIdx); 
  contraIdx                            = laserON & ((trialType == analysisParams.leftCode & locIdx == lidx(2)) | ...
                                                    (trialType == analysisParams.rightCode & locIdx == lidx(1)));
  varPower.vgat.percCorrect_contra(iP) = 100 * sum(trialType(contraIdx) == choice(contraIdx)) / sum(contraIdx); 
  
  laserIdx                             = laserON;
  varPower.vgat.percCorrect_both(iP)   = 100 * sum(trialType(laserIdx) == choice(laserIdx)) / sum(laserIdx); 
  
  [~,varPower.vgat.percCorrect_ctrl_binoErr(iP,:)]   ...
                                       = binofit(sum(lg.trialType(~lg.laserON) == lg.choice(~lg.laserON)),sum(~lg.laserON),1-stdInterval);
  [~,varPower.vgat.percCorrect_ipsi_binoErr(iP,:)]   ...
                                       = binofit(sum(trialType(ipsiIdx) == choice(ipsiIdx)),sum(ipsiIdx),1-stdInterval);
  [~,varPower.vgat.percCorrect_contra_binoErr(iP,:)] ...
                                       = binofit(sum(trialType(contraIdx) == choice(contraIdx)),sum(contraIdx),1-stdInterval);
  [~,varPower.vgat.percCorrect_both_binoErr(iP,:)] ...
                                       = binofit(sum(trialType(laserIdx) == choice(laserIdx)),sum(laserIdx),1-stdInterval);
                                     
  varPower.vgat.excessTravel_ctrl(iP)       = mean(100*lg.excessTravel(~lg.laserON));
  varPower.vgat.excessTravel_ctrl_sem(iP)   = std(100*lg.excessTravel(~lg.laserON))/sqrt(sum(~lg.laserON)-1);
  varPower.vgat.excessTravel_ipsi(iP)       = mean(100*excessTravel(~laserON));
  varPower.vgat.excessTravel_ipsi_sem(iP)   = std(100*excessTravel(ipsiIdx))/sqrt(sum(ipsiIdx)-1);
  varPower.vgat.excessTravel_contra(iP)     = mean(100*excessTravel(contraIdx));
  varPower.vgat.excessTravel_contra_sem(iP) = std(100*excessTravel(contraIdx))/sqrt(sum(contraIdx)-1);
  
  [varPower.stats.vgat.p_percCorrect_ipsi(iP), varPower.stats.vgat.test_percCorrect_ipsi]       ...
                                            = calculatePValue(trialType(ipsiIdx) == choice(ipsiIdx),lg.trialType(~lg.laserON) == lg.choice(~lg.laserON),'boot',0.05);
  [varPower.stats.vgat.p_percCorrect_contra(iP), varPower.stats.vgat.test_percCorrect_contra]   ...
                                            = calculatePValue(trialType(contraIdx) == choice(contraIdx),lg.trialType(~lg.laserON) == lg.choice(~lg.laserON),'boot',0.05);
  [varPower.stats.vgat.p_percCorrect_both(iP), varPower.stats.vgat.test_percCorrect_both, varPower.stats.vgat.sdDiff_percCorrect_both(iP)]   ...
                                            = calculatePValue(trialType(laserON) == choice(laserON),lg.trialType(~lg.laserON) == lg.choice(~lg.laserON),'boot',0.05);
  [varPower.stats.vgat.p_excessTravel_ipsi(iP), varPower.stats.vgat.test_excessTravel_ipsi]     ...
                                            = calculatePValue(excessTravel(ipsiIdx),lg.excessTravel(~lg.laserON),'classic',0.05);
  [varPower.stats.vgat.p_excessTravel_contra(iP), varPower.stats.vgat.test_excessTravel_contra] ...
                                            = calculatePValue(excessTravel(contraIdx),lg.excessTravel(~lg.laserON),'classic',0.05);

  %% ctrl
  choice       = lg_wt.choice(lg_wt.power == varPower.powers(iP));
  trialType    = lg_wt.trialType(lg_wt.power == varPower.powers(iP));
  laserON      = lg_wt.laserON(lg_wt.power == varPower.powers(iP));
  excessTravel = lg_wt.excessTravel(lg_wt.power == varPower.powers(iP));
  locIdx       = cell2mat(lg_wt.galvoPosIdxMaster(lg_wt.power == varPower.powers(iP)));
  
  varPower.stats.ctrl.nTrials_ctrl(iP) = sum(~lg_wt.laserON);
  varPower.stats.ctrl.nTrials_lsr(iP)  = sum(laserON);
  
  varPower.ctrl.percCorrect_ctrl(iP)   = 100 * sum(lg_wt.trialType(~lg_wt.laserON) == lg_wt.choice(~lg_wt.laserON))./sum(~lg_wt.laserON);
  ipsiIdx                              = laserON & ((trialType == analysisParams.leftCode & locIdx == lidx(1)) | ...
                                                    (trialType == analysisParams.rightCode & locIdx == lidx(2)));
  varPower.ctrl.percCorrect_ipsi(iP)   = 100 * sum(trialType(ipsiIdx) == choice(ipsiIdx)) / sum(ipsiIdx); 
  contraIdx                            = laserON & ((trialType == analysisParams.leftCode & locIdx == lidx(2)) | ...
                                                    (trialType == analysisParams.rightCode & locIdx == lidx(1)));
  varPower.ctrl.percCorrect_contra(iP) = 100 * sum(trialType(contraIdx) == choice(contraIdx)) / sum(contraIdx); 
  
  [~,varPower.ctrl.percCorrect_ctrl_binoErr(iP,:)]   ...
                                       = binofit(sum(lg_wt.trialType(~lg_wt.laserON) == lg_wt.choice(~lg_wt.laserON)),sum(~lg_wt.laserON),1-stdInterval);
  [~,varPower.ctrl.percCorrect_ipsi_binoErr(iP,:)]   ...
                                       = binofit(sum(trialType(ipsiIdx) == choice(ipsiIdx)),sum(ipsiIdx),1-stdInterval);
  [~,varPower.ctrl.percCorrect_contra_binoErr(iP,:)] ...
                                       = binofit(sum(trialType(contraIdx) == choice(contraIdx)),sum(contraIdx),1-stdInterval);
                                     
  varPower.ctrl.excessTravel_ctrl(iP)       = mean(100*lg_wt.excessTravel(~lg_wt.laserON));
  varPower.ctrl.excessTravel_ctrl_sem(iP)   = std(100*lg_wt.excessTravel(~lg_wt.laserON))/sqrt(sum(~lg_wt.laserON)-1);
  varPower.ctrl.excessTravel_ipsi(iP)       = mean(100*excessTravel(~laserON));
  varPower.ctrl.excessTravel_ipsi_sem(iP)   = std(100*excessTravel(ipsiIdx))/sqrt(sum(ipsiIdx)-1);
  varPower.ctrl.excessTravel_contra(iP)     = mean(100*excessTravel(contraIdx));
  varPower.ctrl.excessTravel_contra_sem(iP) = std(100*excessTravel(contraIdx))/sqrt(sum(contraIdx)-1);
  
  [varPower.stats.ctrl.p_percCorrect_ipsi(iP), varPower.stats.ctrl.test_percCorrect_ipsi]       ...
                                            = calculatePValue(trialType(ipsiIdx) == choice(ipsiIdx),lg_wt.trialType(~lg_wt.laserON) == lg_wt.choice(~lg_wt.laserON),'boot',0.05);
  [varPower.stats.ctrl.p_percCorrect_contra(iP), varPower.stats.ctrl.test_percCorrect_contra]   ...
                                            = calculatePValue(trialType(contraIdx) == choice(contraIdx),lg_wt.trialType(~lg_wt.laserON) == lg_wt.choice(~lg_wt.laserON),'boot',0.05);
  [varPower.stats.ctrl.p_excessTravel_ipsi(iP), varPower.stats.ctrl.test_excessTravel_ipsi]     ...
                                            = calculatePValue(excessTravel(ipsiIdx),lg_wt.excessTravel(~lg_wt.laserON),'classic',0.05);
  [varPower.stats.ctrl.p_excessTravel_contra(iP), varPower.stats.ctrl.test_excessTravel_contra] ...
                                            = calculatePValue(excessTravel(contraIdx),lg_wt.excessTravel(~lg_wt.laserON),'classic',0.05);                                       
end

varPower.stats.vgat.nTrials_ctrl_mean = mean(varPower.stats.vgat.nTrials_ctrl);
varPower.stats.vgat.nTrials_ctrl_sem  = std(varPower.stats.vgat.nTrials_ctrl)/sqrt(numel(varPower.powers)-1);
varPower.stats.vgat.nTrials_lsr_mean  = mean(varPower.stats.vgat.nTrials_lsr);
varPower.stats.vgat.nTrials_lsr_sem   = std(varPower.stats.vgat.nTrials_lsr)/sqrt(numel(varPower.powers)-1);

varPower.stats.ctrl.nTrials_ctrl_mean = mean(varPower.stats.ctrl.nTrials_ctrl);
varPower.stats.ctrl.nTrials_ctrl_sem  = std(varPower.stats.ctrl.nTrials_ctrl)/sqrt(numel(varPower.powers)-1);
varPower.stats.ctrl.nTrials_lsr_mean  = mean(varPower.stats.ctrl.nTrials_lsr);
varPower.stats.ctrl.nTrials_lsr_sem   = std(varPower.stats.ctrl.nTrials_lsr)/sqrt(numel(varPower.powers)-1);

if spockFlag
  fn = '/jukebox/braininit/Analysis/laserGalvo/V1varPower.mat';
else
  fn = '/Volumes/braininit/Analysis/laserGalvo/V1varPower.mat';  
end
save(fn,'varPower')
