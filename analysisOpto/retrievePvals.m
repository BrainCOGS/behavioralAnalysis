function [pvals,alpha_correct] = retrievePvals(lsrPerf,plotWhat,ptype,alpha_correct)

% pvals = retrievePvals(lsrPerf,plotWhat,ptype)
% builds structure for joint plotting of all locations, returning
% coordinates, effect size and pvalue
% lsrPerf is output structure of analyzeConcatLsrLog()
% plotWhat is one of the following: {'bias','percCorrect','percFinish','speed','trialDur'};
% ptype is 'boot' or 'mouse'
% alpha_correct is the post-FDR significance level
% RETURNS
% pvals.coord [ML AP]
% pvals.p
% pvals.effectSize ([-1 1])

%% defualts
if nargin < 3; ptype = 'boot'; end % 'boot' or 'mouse'
if nargin < 4; alpha_correct = []; end % 'boot' or 'mouse'

%% retrive alpha if necessary
if isempty(alpha_correct)
  switch ptype
    case 'mouse'
      if isfield(lsrPerf{1}.stats.xMouse,['alpha_correct_' plotWhat])
        alpha_correct = eval(['lsrPerf{1}.stats.xMouse.alpha_correct_' plotWhat]);
      else
        alpha_correct = FDRcorrect(lsrPerf,analysisParams.bootalpha,plotWhat);
      end
    case 'boot'
      if isfield(lsrPerf{1}.FDR,['alpha_correct_' plotWhat])
        alpha_correct = eval(['lsrPerf{1}.FDR.alpha_correct_' plotWhat]);
      else
        if ~isempty(strfind(plotWhat,'viewDecode'))
          alpha_correct = lsrPerf{1}.FDR.alpha_correct_viewDecode;
        else
          alpha_correct = analysisParams.bootalpha;%FDRcorrect(lsrPerf,analysisParams.bootalpha,plotWhat);
        end
      end
  end
end

%% retrieve effect size  and significance
if iscell(lsrPerf{1}.info.grid) % grids w/simulatenous locations
  
  pvals.coord      = [];
  pvals.p          = [];
  pvals.effectSize = [];
  for ll = 1:length(lsrPerf)
    if ll == 1; nloc = size(lsrPerf{1}.info.grid{ll},1); end
    pvals.coord      = [pvals.coord; lsrPerf{1}.info.grid{ll}];
    [p,es]           = getPlocation(lsrPerf,ll,plotWhat,ptype);
    pvals.p          = [pvals.p; repmat(p,[nloc 1])];
    pvals.effectSize = [pvals.effectSize; repmat(es,[nloc 1])];
  end
else
  [~,mgr] = getMasterLocationIdx([]);
  for ll = 1:length(lsrPerf)
    pvals.coord(ll,:)                         = mgr(lsrPerf{ll}.lidxMaster,:);
    [pvals.p(ll),pvals.effectSize(ll)]        = getPlocation(lsrPerf,ll,plotWhat,ptype);
    pvals.effectSize(isnan(pvals.effectSize)) = 0;
    pvals.p(isnan(pvals.p))                   = 1;
  end
end

end

function [p,es] = getPlocation(lsrPerf,ll,plotWhat,ptype)

switch plotWhat
  case {'bias','bias_abs','percCorrect','percCorrect_R','percCorrect_L', ...
      'percCorrect_easy','percCorrect_hard','lapse','percFinish'}
    switch ptype
      case 'boot'
        es = .01*(lsrPerf{ll}.lsr.(plotWhat)-lsrPerf{ll}.ctrl.(plotWhat));
        p  = lsrPerf{ll}.stats.bootPerf.(['p_' plotWhat]);
        
      case 'mouse'
        es = .01*(lsrPerf{ll}.lsr.(plotWhat)-lsrPerf{ll}.ctrl.(plotWhat));
        p  = lsrPerf{ll}.stats.xMouse.(['p_' plotWhat]);
    end
    
  case 'percCorrect_hard_norm'
    nt_ctrl = lsrPerf{ll}.stats.ctrlTrialsPerMouse;
    nt_lsr  = lsrPerf{ll}.stats.lsrTrialsPerMouse;
    ctrl    = lsrPerf{ll}.ctrl.percCorrect_hard_mouse;
    lsr     = lsrPerf{ll}.lsr.percCorrect_hard_mouse;
    nanidx  = isnan(ctrl) | isnan(lsr);
    ctrl(nanidx)    = [];
    lsr(nanidx)     = [];
    nt_ctrl(nanidx) = [];
    nt_lsr(nanidx)  = [];
    es      = .01 * (((lsr*nt_lsr')./sum(nt_lsr))-((ctrl*nt_ctrl')./sum(nt_ctrl)) ./ ((ctrl*nt_ctrl')./sum(nt_ctrl)-.5));
    p  = lsrPerf{ll}.stats.bootPerf.p_percCorrect_hard;
    
  case 'percCorrect_easy_norm'
    nt_ctrl = lsrPerf{ll}.stats.ctrlTrialsPerMouse;
    nt_lsr  = lsrPerf{ll}.stats.lsrTrialsPerMouse;
    ctrl    = lsrPerf{ll}.ctrl.percCorrect_easy_mouse;
    lsr     = lsrPerf{ll}.lsr.percCorrect_easy_mouse;
    nanidx  = isnan(ctrl) | isnan(lsr);
    ctrl(nanidx)    = [];
    lsr(nanidx)     = [];
    nt_ctrl(nanidx) = [];
    nt_lsr(nanidx)  = [];
    es      = .01 * (((lsr*nt_lsr')./sum(nt_lsr))-((ctrl*nt_ctrl')./sum(nt_ctrl)) ./ ((ctrl*nt_ctrl')./sum(nt_ctrl)-.5));
    p  = lsrPerf{ll}.stats.bootPerf.p_percCorrect_easy;
    
  case 'percCorrect_norm'
    switch ptype
      case 'boot'
        es = .1*(lsrPerf{ll}.lsr.(plotWhat));
        p  = lsrPerf{ll}.stats.bootPerf.(['p_' plotWhat]);
        
      case 'mouse'
        es = .01*mean(lsrPerf{ll}.lsr.([plotWhat '_mouse']));
        p  = lsrPerf{ll}.stats.xMouse.(['p_' plotWhat]);
    end
    
  case {'speed','dur'}
    switch ptype
      case 'boot'
        p  = lsrPerf{ll}.stats.(['p_' plotWhat]);
        es = (lsrPerf{ll}.lsr.(plotWhat)-lsrPerf{ll}.ctrl.(plotWhat))/lsrPerf{ll}.ctrl.(plotWhat);
        
      case 'mouse'
        p  = lsrPerf{ll}.stats.xMouse.(['p_' plotWhat]);
        es = mean((lsrPerf{ll}.lsr.([plotWhat '_mouse'])-lsrPerf{ll}.ctrl.([plotWhat '_mouse']))./lsrPerf{ll}.ctrl.([plotWhat '_mouse']));
    end
    
  case 'excessTravel'
    switch ptype
      case 'boot'
        p  = lsrPerf{ll}.stats.(['p_' plotWhat]);
        es = .01.*(lsrPerf{ll}.lsr.(plotWhat)-lsrPerf{ll}.ctrl.(plotWhat));
        
      case 'mouse'
        p  = lsrPerf{ll}.stats.xMouse.(['p_' plotWhat]);
        es = .01.*mean(lsrPerf{ll}.lsr.([plotWhat '_mouse'])-lsrPerf{ll}.ctrl.([plotWhat '_mouse']));
    end
    
  case 'unusualEvents'
    switch ptype
      case 'boot'
        p  = lsrPerf{ll}.stats.bootMotor.(['p_' plotWhat]);
        es = .01*(lsrPerf{ll}.lsr.(plotWhat)-lsrPerf{ll}.ctrl.(plotWhat))/lsrPerf{ll}.ctrl.(plotWhat);
        
      case 'mouse'
        p  = lsrPerf{ll}.stats.xMouse.(['p_' plotWhat]);
        es = .01*mean((lsrPerf{ll}.lsr.([plotWhat '_mouse'])-lsrPerf{ll}.ctrl.([plotWhat '_mouse']))./lsrPerf{ll}.ctrl.([plotWhat '_mouse']));
    end
    
  case 'slope'
    p  = lsrPerf{ll}.stats.bootPerf.p_slope;
    es = .001*((100*lsrPerf{ll}.lsr.psychometric.fit.slope-100*lsrPerf{ll}.ctrl.psychometric.fit.slope)/(100*lsrPerf{ll}.ctrl.psychometric.fit.slope));
    
  case 'viewAngSD'
    p  = lsrPerf{ll}.viewAngle.p_viewAngSD;
    es = (lsrPerf{ll}.viewAngle.viewAngSD_lsr-lsrPerf{ll}.viewAngle.viewAngSD_ctrl)/lsrPerf{ll}.viewAngle.viewAngSD_ctrl;
    
  case 'viewDecode'
    p  = lsrPerf{ll}.viewAngle.p_decodeAcc_overall;
    es = mean(lsrPerf{ll}.viewAngle.decodeAcc_diff);
    
  case 'viewDecodeCueHalf1'
    nb = numel(lsrPerf{ll}.viewAngle.decodeAcc_diff);
    p  = lsrPerf{ll}.viewAngle.p_decodeAcc_cueHalf1;
    es = mean(lsrPerf{ll}.viewAngle.decodeAcc_diff(1:round(nb/3)));
    
  case 'viewDecodeCueHalf2'
    nb = numel(lsrPerf{ll}.viewAngle.decodeAcc_diff);
    p  = lsrPerf{ll}.viewAngle.p_decodeAcc_cueHalf2;
    es = mean(lsrPerf{ll}.viewAngle.decodeAcc_diff(round(nb/3):round(nb/3)*2));
    
  case 'viewDecodeCue'
    nb = numel(lsrPerf{ll}.viewAngle.decodeAcc_diff);
    p  = lsrPerf{ll}.viewAngle.p_decodeAcc_cue;
    es = mean(lsrPerf{ll}.viewAngle.decodeAcc_diff(1:round(nb/3)*2));
    
  case 'viewDecodeMem'
    nb = numel(lsrPerf{ll}.viewAngle.decodeAcc_diff);
    p  = lsrPerf{ll}.viewAngle.p_decodeAcc_mem;
    es = mean(lsrPerf{ll}.viewAngle.decodeAcc_diff(round(nb/3)*2:end));
    
  case 'logRegSlope'
    p  = lsrPerf{ll}.lsr.logisticReg.slope_p;
    es = lsrPerf{ll}.lsr.logisticReg.slope_diff_mean;
    
  case 'logRegDecayIndex'
    p  = lsrPerf{ll}.lsr.logisticReg.decayIndex_p;
    es = lsrPerf{ll}.lsr.logisticReg.decayIndex_diff_mean;
    
  case 'viewOverlap'
    p  = lsrPerf{ll}.stats.bootPerf.p_viewOverlap;
    es = nanmean(lsrPerf{ll}.lsr.trajectory.viewOverlap_coarseBinned(1:end-1)-lsrPerf{ll}.ctrl.trajectory.viewOverlap_coarseBinned(1:end-1));
end
end