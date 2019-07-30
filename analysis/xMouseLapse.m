function lapse = xMouseLapse(lg)

% lapse = xMouseLapse(lg)
% average per-mouse lapse rate
% lg is flattened behavior log

mice        = unique(lg.mouseID);
lapse       = nan(numel(mice), 1);

for iMouse = 1:numel(mice)
  [choice, trialtype, delta]  = selectMouseTrials(lg, mice(iMouse), 'choice', 'trialType', 'nCues_RminusL');
  lapse(iMouse)               = 100* sum(choice(abs(delta)>=10) ~= trialtype(abs(delta)>=10))/sum(abs(delta)>=10);
end