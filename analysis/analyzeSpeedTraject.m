function sta = analyzeSpeedTraject(lgfile,cfg)

% sta = analyzeSpeedTraject(lgfile,cfg)
% analyzes speed and trajectories during vr navigation
% lg file can be either the flattened log or its path
%
% Lucas Pinto (lpinto@princeton.edu)


%% load if lgfile is a path
if ischar(lgfile)
  load(lgfile,'lg')
else
  lg = lgfile;
end

if nargin < 2
  cfg = struct([]);
end

%% default analysis params
cfg   = populateCfg(cfg);

%% just in case it hasn't been done before
lg    = cleanupConcatLog(lg, cfg.minNumTrials);

%% loop through mice
mice                              = unique(lg.mouseID);

sta.speedAvg_mouse                = nan(numel(mice),1);
sta.perfAvg_mouse                 = nan(numel(mice),1);
sta.speedStd_withinSess_mouse     = nan(numel(mice),1);
sta.speedVsPerf_sess_mouse_corr   = nan(numel(mice),1);
sta.speedVsPerf_sess_mouse_pval   = nan(numel(mice),1);
sta.speedStd_xSess_mouse          = nan(numel(mice),1);
sta.viewAngleStd_withinSess_mouse = nan(numel(mice),1);
sta.viewAngleStd_xSess_mouse      = nan(numel(mice),1);


for iMouse = 1:numel(mice)
  [speedStem, sessionID, pos, choice, trialType]    ...
                             = selectMouseTrials(lg, mice(iMouse), 'speedStem', 'sessionID', 'pos', 'choice', 'trialType');
                      
  sta.speedAvg_mouse(iMouse) = nanmean(speedStem);
  sta.perfAvg_mouse(iMouse)  = sum(choice == trialType)/numel(choice) * 100;
  
  sess                       = unique(sessionID);
  
  speedSess_avg              = nan(numel(sess),1);
  speedSess_std              = nan(numel(sess),1);
  perfSess                   = nan(numel(sess),1);
  viewSess_avg               = nan(numel(cfg.posBins),numel(sess),2);
  viewSess_std               = nan(numel(sess),1);
  
  for iSess = 1:numel(sess)
    iSpeed  = speedStem(sessionID == sess(iSess));
    iPos    = pos      (sessionID == sess(iSess));
    iChoice = choice   (sessionID == sess(iSess));
    iTtype  = trialType(sessionID == sess(iSess));
    
    % speed
    if numel(iSpeed) < cfg.minNumTrialSess; continue; end
    
    speedSess_avg(iSess)     = nanmean(iSpeed);
    speedSess_std(iSess)     = nanstd(iSpeed);   
    perfSess(iSess)          = sum(iChoice == iTtype)/numel(iChoice)*100;   
    
    % do angle by Y pos, taking into account the first time that y is
    % reached
    viewAngle                = sampleViewAngleVsY(iPos, cfg.posBins);
    viewAngleR               = viewAngle(:,iChoice==analysisParams.rightCode);
    viewAngleL               = viewAngle(:,iChoice==analysisParams.leftCode);
    viewSess_avg(:,iSess,1)  = nanmean(viewAngleR,2);
    viewSess_avg(:,iSess,2)  = nanmean(viewAngleL,2);
    viewSess_std(iSess)      = nanmean(nanmean([nanstd(viewAngleR,0,2) nanstd(viewAngleL,0,2)],2));
    
  end
  
  sta.speedStd_withinSess_mouse(iMouse)     = nanmean(speedSess_std);
  sta.speedStd_xSess_mouse(iMouse)          = nanstd(speedSess_avg);
  sta.viewAngleStd_withinSess_mouse(iMouse) = nanmean(viewSess_std);
  sta.viewAngleStd_xSess_mouse(iMouse)      = nanmean(nanmean([nanstd(viewSess_avg(:,:,1),0,2) nanstd(viewSess_avg(:,:,2),0,2)],2));

  % speed vs perf correlation on sessions
  [sta.speedVsPerf_sess_mouse_corr(iMouse), sta.speedVsPerf_sess_mouse_pval(iMouse)]                        ...
                                            = corr(speedSess_avg(~isnan(speedSess_avg) & ~isnan(perfSess)), ...
                                                   perfSess     (~isnan(speedSess_avg) & ~isnan(perfSess)));
end

%% collect stats
sta.stats.ngoodmice                   = sum(~isnan(sta.speedStd_withinSess_mouse));
nnorm                                 = sqrt(sta.stats.ngoodmice-1);
sta.stats.speed_avg                   = nanmean(sta.speedAvg_mouse);
sta.stats.speed_sem                   = nanstd(sta.speedAvg_mouse)/nnorm;
sta.stats.speedStd_withinSess_avg     = nanmean(sta.speedStd_withinSess_mouse);
sta.stats.speedStd_withinSess_sem     = nanstd(sta.speedStd_withinSess_mouse)/nnorm;
sta.stats.speedStd_xSess_avg          = nanmean(sta.speedStd_xSess_mouse);
sta.stats.speedStd_xSess_sem          = nanstd(sta.speedStd_xSess_mouse)/nnorm;

[sta.stats.speedVsPerf_corr,sta.stats.speedVsPerf_pval] ...
                                      = corr(sta.speedAvg_mouse, sta.perfAvg_mouse);

sta.stats.viewAngleStd_withinSess_avg = nanmean(sta.viewAngleStd_withinSess_mouse);
sta.stats.viewAngleStd_withinSess_sem = nanstd(sta.viewAngleStd_withinSess_mouse)/nnorm;
sta.stats.viewAngleStd_xSess_avg      = nanmean(sta.viewAngleStd_xSess_mouse);
sta.stats.viewAngleStd_xSess_sem      = nanstd(sta.viewAngleStd_xSess_mouse)/nnorm;

end

%% cfg defaults
function cfg = populateCfg(cfg)

if ~isfield(cfg,'minNumTrials')
  cfg(1).minNumTrials      = 1000;
end
if ~isfield(cfg,'minNumTrialSess')
  cfg(1).minNumTrialSess   = 0;
end
if ~isfield(cfg,'posBins')
  cfg(1).posBins           = 0:300;
end

end