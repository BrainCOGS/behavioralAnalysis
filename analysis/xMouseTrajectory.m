function traject = xMouseTrajectory(log)

% traject = xMouseTrajectory(log)
% average per-mouse view angle trajectories
% lg is flattened behavior log

% delete aborted / incomplete trials
idx = find(log.choice~=analysisParams.rightCode & log.choice~=analysisParams.leftCode);
log.choice(idx)         = [];
log.pos(idx)            = [];
log.excessTravel(idx)   = [];
log.sessionID(idx)      = [];
log.cuePos_R(idx)       = [];
log.cuePos_L(idx)       = [];
log.mouseID(idx)        = [];
log.speedStem(idx)      = [];

% overall first
traject        = analyzeTrajectories(log.pos,log.choice,log.excessTravel);
traject.vaVSev = viewAngleVSevidence(log.pos,log.cuePos_R,log.cuePos_L,log.choice,log.excessTravel);

traject.mice   = unique(log.mouseID);
traject.nmice  = length(traject.mice);


tempViewR = []; tempViewL = []; tempOverlap = [];
traject.speedAvg = []; traject.speedSD_xTrial = []; traject.speedSD_xSess = [];
traject.corr_instVaVSev = []; traject.corr_cueHalf1VaVSev = []; traject.corr_cueAllVaVSev = [];
traject.corr_cueAllVaVSevL = []; traject.corr_cueAllVaVSevR = [];
tempinstVaVSev = [];  tempch1VaVSev = []; tempch2VaVSev = []; tempcaVaVSev = [];
tempch1VaVSevR = []; tempch2VaVSevR = []; tempcaVaVSevR = [];
tempch1VaVSevL = []; tempch2VaVSevL = []; tempcaVaVSevL = [];
for ii = 1:traject.nmice
  traject.mouse(ii).nTrials   = sum(log.mouseID == traject.mice(ii));
  traject.mouse(ii).traject   = analyzeTrajectories(log.pos(log.mouseID == traject.mice(ii)),log.choice(log.mouseID == traject.mice(ii)),log.excessTravel(log.mouseID == traject.mice(ii)));
  traject.mouse(ii).speedAvg  = nanmean(log.speedStem(log.mouseID == traject.mice(ii) & log.excessTravel < .1));
  traject.mouse(ii).speedSD   = nanstd(log.speedStem(log.mouseID == traject.mice(ii) & log.excessTravel < .1));
  traject.mouse(ii).vaVSev    = viewAngleVSevidence(log.pos(log.mouseID == traject.mice(ii)),log.cuePos_R(log.mouseID == traject.mice(ii)),log.cuePos_L(log.mouseID == traject.mice(ii)),log.choice(log.mouseID == traject.mice(ii)),log.excessTravel(log.mouseID == traject.mice(ii)));
  
  % go through individual sessions for speed stereotypy
  sess  = unique(log.sessionID(log.mouseID == traject.mice(ii)));
  avgsp = []; sessSD = [];
  for jj = 1:length(sess)
    avgsp(jj)  = nanmean(log.speedStem(log.mouseID == traject.mice(ii) & log.sessionID == sess(jj) & log.excessTravel < .1));
    sessSD(jj) = nanstd(log.speedStem(log.mouseID == traject.mice(ii) & log.sessionID == sess(jj) & log.excessTravel < .1));
  end
  
  traject.mouse(ii).speedSD_xTrial  = nanmean(sessSD);
  traject.mouse(ii).speedSD_xSess   = nanstd(avgsp);
  
  if traject.mouse(ii).nTrials > 200
    tempViewR(:,:,end+1)                = traject.mouse(ii).traject.viewR_coarseBinned;
    tempViewL(:,:,end+1)                = traject.mouse(ii).traject.viewL_coarseBinned;
    tempOverlap(:,end+1)                = traject.mouse(ii).traject.viewOverlap_coarseBinned;
    traject.speedAvg(end+1)             = traject.mouse(ii).speedAvg;
    traject.speedSD_xTrial(end+1)       = traject.mouse(ii).speedSD_xTrial;
    traject.speedSD_xSess(end+1)        = traject.mouse(ii).speedSD_xSess;
    traject.corr_instVaVSev(end+1)      = corr(traject.mouse(ii).vaVSev.RmL',traject.mouse(ii).vaVSev.inst_angleVSevidence_mean');
    traject.corr_cueHalf1VaVSev(end+1)  = corr(traject.mouse(ii).vaVSev.RmL(~isnan(traject.mouse(ii).vaVSev.cueHalf1_angleVSevidence_mean))',...
      traject.mouse(ii).vaVSev.cueHalf1_angleVSevidence_mean(~isnan(traject.mouse(ii).vaVSev.cueHalf1_angleVSevidence_mean))');
    traject.corr_cueAllVaVSev(end+1)    = corr(traject.mouse(ii).vaVSev.RmL(~isnan(traject.mouse(ii).vaVSev.cueAll_angleVSevidence_mean))',...
      traject.mouse(ii).vaVSev.cueAll_angleVSevidence_mean(~isnan(traject.mouse(ii).vaVSev.cueAll_angleVSevidence_mean))');
    traject.corr_cueAllVaVSevR(end+1)   = corr(abs(traject.mouse(ii).vaVSev.RmL(~isnan(traject.mouse(ii).vaVSev.cueAll_angleVSevidence_R_mean)))',...
      abs(traject.mouse(ii).vaVSev.cueAll_angleVSevidence_R_mean(~isnan(traject.mouse(ii).vaVSev.cueAll_angleVSevidence_R_mean))'));
    traject.corr_cueAllVaVSevL(end+1)   = corr(abs(traject.mouse(ii).vaVSev.RmL(~isnan(traject.mouse(ii).vaVSev.cueAll_angleVSevidence_L_mean)))',...
      abs(traject.mouse(ii).vaVSev.cueAll_angleVSevidence_L_mean(~isnan(traject.mouse(ii).vaVSev.cueAll_angleVSevidence_L_mean))'));
    tempinstVaVSev(:,end+1)             = traject.mouse(ii).vaVSev.inst_angleVSevidence_mean';
    tempch1VaVSev(:,end+1)              = traject.mouse(ii).vaVSev.cueHalf1_angleVSevidence_mean';
    tempch2VaVSev(:,end+1)              = traject.mouse(ii).vaVSev.cueHalf2_angleVSevidence_mean';
    tempcaVaVSev(:,end+1)               = traject.mouse(ii).vaVSev.cueAll_angleVSevidence_mean';
    tempch1VaVSevR(:,end+1)             = traject.mouse(ii).vaVSev.cueHalf1_angleVSevidence_R_mean';
    tempch2VaVSevR(:,end+1)             = traject.mouse(ii).vaVSev.cueHalf2_angleVSevidence_R_mean';
    tempcaVaVSevR(:,end+1)              = traject.mouse(ii).vaVSev.cueAll_angleVSevidence_R_mean';
    tempch1VaVSevL(:,end+1)             = traject.mouse(ii).vaVSev.cueHalf1_angleVSevidence_L_mean';
    tempch2VaVSevL(:,end+1)             = traject.mouse(ii).vaVSev.cueHalf2_angleVSevidence_L_mean';
    tempcaVaVSevL(:,end+1)              = traject.mouse(ii).vaVSev.cueAll_angleVSevidence_L_mean';
  end
end

traject.ViewRmean               = nanmean(tempViewR,3);
traject.ViewLmean               = nanmean(tempViewL,3);
traject.overlapMean             = nanmean(tempOverlap,2);
traject.instAngleVSev_mice      = tempinstVaVSev;
traject.cueHalf1AngleVSev_mice  = tempch1VaVSev;
traject.cueAllAngleVSev_mice    = tempcaVaVSev;
traject.instAngleVSevMean       = nanmean(tempinstVaVSev,2);
traject.cueHalf1AngleVSevMean   = nanmean(tempch1VaVSev,2);
traject.cueHalf2AngleVSevMean   = nanmean(tempch2VaVSev,2);
traject.cueAllAngleVSevMean     = nanmean(tempcaVaVSev,2);
traject.cueHalf1AngleVSevMeanR  = nanmean(tempch1VaVSevR,2);
traject.cueAllAngleVSevMeanR    = nanmean(tempcaVaVSevR,2);
traject.cueHalf1AngleVSevMeanL  = nanmean(tempch1VaVSevL,2);
traject.cueAllAngleVSevMeanL    = nanmean(tempcaVaVSevL,2);
traject.cueHalf2AngleVSevMeanR  = nanmean(tempch2VaVSevR,2);
traject.cueHalf2AngleVSevMeanL  = nanmean(tempch2VaVSevL,2);
traject.ViewRsem                = nanstd(tempViewR,0,3)./sqrt(size(tempViewR,3)-1);
traject.ViewLsem                = nanstd(tempViewL,0,3)./sqrt(size(tempViewR,3)-1);
traject.overlapSem              = nanstd(tempOverlap,0,2)./sqrt(size(tempViewR,3)-1);
traject.instAngleVSevSem        = nanstd(tempinstVaVSev,0,2)./sqrt(size(tempcaVaVSevR,2)-1);
traject.cueHalf1AngleVSevSem    = nanstd(tempch1VaVSev,0,2)./sqrt(size(tempcaVaVSevR,2)-1);
traject.cueHalf2AngleVSevSem    = nanstd(tempch2VaVSev,0,2)./sqrt(size(tempcaVaVSevR,2)-1);
traject.cueAllAngleVSevSem      = nanstd(tempcaVaVSev,0,2)./sqrt(size(tempcaVaVSevR,2)-1);
traject.cueHalf1AngleVSevSemR   = nanstd(tempch1VaVSevR,0,2)./sqrt(size(tempcaVaVSevR,2)-1);
traject.cueHalf2AngleVSevSemR   = nanstd(tempch2VaVSevR,0,2)./sqrt(size(tempcaVaVSevR,2)-1);
traject.cueAllAngleVSevSemR     = nanstd(tempcaVaVSevR,0,2)./sqrt(size(tempcaVaVSevR,2)-1);
traject.cueHalf1AngleVSevSemL   = nanstd(tempch1VaVSevL,0,2)./sqrt(size(tempcaVaVSevR,2)-1);
traject.cueHalf2AngleVSevSemL   = nanstd(tempch2VaVSevL,0,2)./sqrt(size(tempcaVaVSevR,2)-1);
traject.cueAllAngleVSevSemL     = nanstd(tempcaVaVSevL,0,2)./sqrt(size(tempcaVaVSevR,2)-1);