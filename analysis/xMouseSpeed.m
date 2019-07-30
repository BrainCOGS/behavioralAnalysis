function speed = xMouseSpeed(lg,mazeSection)

% speed = xMouseSpeed(lg,mazeSection)
% average per-mouse running speed
% lg is flattened behavior log
% mazeSection is string specifying analysis window: 'all','cue' or 'delay'

if nargin < 2; mazeSection = 'all'; end

mice        = unique(lg.mouseID);
speed       = nan(numel(mice), 1);

for iMouse = 1:numel(mice)
  switch mazeSection
    case 'all'
      speedStem     = selectMouseTrials(lg, mice(iMouse), 'speedStem');
      speed(iMouse) = nanmean(speedStem);
  
    case {'cue','delay'}
      [pos,time]    = selectMouseTrials(lg, mice(iMouse), 'pos', 'time');
      ntrials       = numel(pos);
      speedTrials   = nan(ntrials,1);
      for itrial = 1:ntrials
        switch mazeSection
          case 'cue'
            startidx    = find(pos{itrial}(:,2) >= 0,   1, 'first');
            endidx      = find(pos{itrial}(:,2) <= 200, 1, 'last');
          case 'delay'
            startidx    = find(pos{itrial}(:,2) >  200, 1, 'first');
            endidx      = find(pos{itrial}(:,2) <  300, 1, 'last');
        end
        displ       = diff(pos{itrial}(startidx:endidx,1:2));
        displ       = sum(sqrt(sum(displ.^2,2)));
        elt         = time{itrial}(endidx) - time{itrial}(startidx);
        speedTrials(itrial) = displ / elt;
      end
      speed(iMouse) = nanmean(speedTrials);
      
  end
end