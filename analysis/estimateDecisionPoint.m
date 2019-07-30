function [decisionPoint,net_ev,viewAngle,deriv,trialIdx,cfg] = estimateDecisionPoint(lg,cfg,viewAngle,selectTrials,method)

% [decisionPoint,net_ev,viewAngle,deriv,trialIdx] = estimateDecisionPoint(lg,cfg,viewAngle,selectTrials,method)

%% analysis config
if nargin < 2 || isempty(cfg)
  cfg.posBins   = 0:300;
  cfg.basel     = 75;%50
  cfg.nSD       = 3;%%1.96;
  cfg.nConsec   = 4;
  cfg.trialType = 'all';
  cfg.secDeriv  = false;
  cfg.velTh     = 0.004; 
  cfg.smoothWin = 12;
end
if nargin < 3 || isempty(viewAngle);    viewAngle    = [];   end
if nargin < 4 || isempty(selectTrials); selectTrials = true; end
if nargin < 5 || isempty(method);       method = 'velocity'; end % velocity or derivative

cfg.method = method;

%% bin viewAngles for desired trial type
if isempty(viewAngle) && selectTrials
  trialIdx    = getTrialIdx(lg,cfg.trialType);
  viewAngle   = sampleViewAngleVsY(lg.pos(trialIdx), cfg.posBins);
elseif isempty(viewAngle) && ~selectTrials
  trialIdx    = 1:numel(lg.choice);
  viewAngle   = sampleViewAngleVsY(lg.pos(trialIdx), cfg.posBins);
else
  trialIdx    = 1:size(viewAngle,2);
  cfg.posBins = linspace(0,300,size(viewAngle,1));
end

%% for each trial, estimate decision point
ntrials       = numel(trialIdx);
decisionPoint = nan(ntrials,1);
deriv         = nan(numel(cfg.posBins),ntrials);
net_ev        = decisionPoint;
choice        = -(lg.choice .* 2 - 1);
for iTrial = 1:ntrials
  
  switch cfg.method
    
    case 'derivative'
      % take derivative and find position at which dtheta/dy goes beyond cfg.nSD
      % standard deviations from the mean during the baseline period cfg.basel,
      % for cfg.nConsec spatial bins
      deriv(:,iTrial) = [0; diff(viewAngle(:,iTrial))]; 
      if cfg.secDeriv; deriv(:,iTrial) = [0; diff(deriv(:,iTrial))]; end
      basel           = 1:find(cfg.posBins<=cfg.basel,1,'last');
      avg             = mean(deriv(basel,iTrial));
      sd              = std(deriv(basel,iTrial));
      thd1            = avg + cfg.nSD*sd;
      thd2            = avg - cfg.nSD*sd;
      if lg.choice(trialIdx(iTrial)) == analysisParams.rightCode
        islarger      = deriv(:,iTrial) > max([thd1 thd2]);
      else
        islarger      = deriv(:,iTrial) < min([thd1 thd2]);
      end
      for iConsec = 1:cfg.nConsec
        islarger(:,end+1) = [islarger(2:end,end); false];
      end
      dP              = cfg.posBins(find(sum(islarger,2) == cfg.nConsec,1,'first'));
      
    case 'velocity'
      % per Harvey et al 2012, consider turn to be last point with x axis
      % velocity of less than cfg.velTh rotations/sec
      switch analysisParams.genotypes{lg.genotypeID(iTrial)} % rig calibartion
        case {'vgat','wt'}
          cal = analysisParams.dotsPerRev;
        otherwise
          cal = widefieldParams.dotsPerRev;
      end
      
      dX    = getXvelocity(lg.sensorDots{iTrial}(1:size(lg.pos{iTrial},1),:),cal);
      dXabs = smooth(abs(dX),cfg.smoothWin);
      dPidx = find(dXabs < cfg.velTh & lg.pos{iTrial}(:,2) > cfg.basel & lg.pos{iTrial}(:,2) < 295, 1, 'last');
      dP    = lg.pos{iTrial}(dPidx,2);
      
  end
  
  if ~isempty(dP); decisionPoint(iTrial) = dP; end
  
  % absolute amount of evidence at decision point
  nR                    = sum(lg.cuePos_R{trialIdx(iTrial)}-10 <= decisionPoint(iTrial));
  nL                    = sum(lg.cuePos_L{trialIdx(iTrial)}-10 <= decisionPoint(iTrial));
  net_ev(iTrial)        = abs(nR-nL);
  
end

end

%% rotation from sensor data
function dX = getXvelocity(sensorDots,cal)

dX = sensorDots(:,4);
dT = sensorDots(:,5)./1000;
dX = 1/cal .* dX ./ dT;

end
