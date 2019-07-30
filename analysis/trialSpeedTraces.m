function [trialSpeed, timeAxis, epochVec] = trialSpeedTraces(pos,time,cfg)

% [trialSpeed, timeAxis, epochVec] = trialSpeedTraces(pos,time,cfg)

%% default analysis params
if nargin < 3
  cfg = struct([]);
end

cfg   = populateCfg(cfg);

%% for each trial, build a vector with sliding average of speed and equivalent time points
% align things to cue period onset

ntrials    = numel(pos);
trialSpeed = cell(1,ntrials);
timeAxis   = cell(1,ntrials);
epochVec   = cell(1,ntrials);

for iTrial = 1:ntrials
  x     = pos{iTrial}(:,1);
  y     = pos{iTrial}(:,2);
  t     = time{iTrial}(1:numel(x));
  zerot = t(find(y>=0,1,'first'));
  t     = t - zerot; % time = 0 is beggining of cue period
  
  tbins = min(t):cfg.wsize:max(t);
  speed = nan(1,numel(tbins)-1);
  epoch = speed;
  
  % calculate speed in each bin
  for iBin = 1:numel(tbins)-1
    idx         = find(t >= tbins(iBin) & t < tbins(iBin+1));
    if isempty(idx); continue; end
    yd          = y(idx(end)) - y(idx(1)); % y displacement 
    xd          = x(idx(end)) - x(idx(1)); % x displacement
    speed(iBin) = sqrt(yd^2+xd^2) / cfg.wsize; % euclidian displacement / time step
    
    ym          = nanmean(y(idx));
    
    if ym < 0
      epoch(iBin) = 1;
    elseif ym >= 0 && ym <= 200
      epoch(iBin) = 2;
    elseif ym > 200 && ym < 280
      epoch(iBin) = 3;
    else
      epoch(iBin) = 4;
    end
    
  end
  
  trialSpeed{iTrial} = speed;
  epochVec{iTrial}   = epoch;
  timeAxis{iTrial}   = toBinCenters(tbins);
end


end
%% cfg defaults
function cfg = populateCfg(cfg)

if ~isfield(cfg,'lStart')
  cfg(1).lStart      = -30;
else
  if cfg.lStart > 0; cfg.lStart = - cfg.lStart; end
end
if ~isfield(cfg,'lCue')
  cfg(1).lCue        = 200;
end
if ~isfield(cfg,'lMemory')
  cfg(1).lMemory     = 100;
end
if ~isfield(cfg,'lMaze')
  cfg(1).lMaze       = cfg.lCue + cfg.lMemory;
end
if ~isfield(cfg,'turnTh')
  cfg(1).turnTh      = 0.5;
end
if ~isfield(cfg,'turnThConsecPoints')
  cfg(1).turnThConsecPoints      = 10;
end
if ~isfield(cfg,'wsize')
  cfg(1).wsize       = .25;
end

end