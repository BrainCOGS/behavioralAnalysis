function turnOnset = xMouseTurnOnset(lg,cfg)

% turnOnset = xMouseTurnOnset(lg,cfg)
% average per-mouse turn onset point
% lg is flattened behavior log
% cfg is optional analysis parameter structure

if nargin < 2; cfg = []; end
if isempty(cfg)
  cfg.posBins   = 0:300;
  cfg.basel     = 75;%50
  cfg.nSD       = 3;%%1.96;
  cfg.nConsec   = 4;
  cfg.trialType = 'all';
  cfg.secDeriv  = false;
  cfg.velTh     = 0.004; 
  cfg.smoothWin = 12;
  cfg.method    = 'derivative';
end

mice        = unique(lg.mouseID);
turnOnset   = nan(numel(mice), 1);

for iMouse = 1:numel(mice)
  sublg             = selectMouseTrials(lg, mice(iMouse));
  thisturn          = estimateDecisionPoint(sublg,cfg,[],false,cfg.method);
  turnOnset(iMouse) = nanmean(thisturn);
end