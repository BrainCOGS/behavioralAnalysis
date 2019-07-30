function accel = xMouseAcceleration(lg,mazeSection)

% accel = xMouseAcceleration(lg,mazeSection)
% average per-mouse acceleration 
% lg is flattened behavior log
% mazeSection is range for calculation in cm, default [250 300]

if nargin < 2; mazeSection = [250 300]; end

mice        = unique(lg.mouseID);
accel       = nan(numel(mice), 1);

for iMouse = 1:numel(mice)
  [pos,time]    = selectMouseTrials(lg, mice(iMouse), 'pos', 'time');
  ntrials       = numel(pos);
  accelTrials   = nan(ntrials,1);
  
  for itrial = 1:ntrials
    startidx    = find(pos{itrial}(:,2) >= mazeSection(1),   1, 'first');
    endidx      = find(pos{itrial}(:,2) <  mazeSection(2), 1, 'last');
    displ       = diff(pos{itrial}(startidx:endidx,1:2));
    displ       = sqrt(sum(displ.^2,2));
    elt         = time{itrial}(startidx:endidx);
    speed       = displ ./ diff(elt);
    accelTrials(itrial) = mean(diff(speed));
  end
  accel(iMouse) = nanmean(accelTrials);

end