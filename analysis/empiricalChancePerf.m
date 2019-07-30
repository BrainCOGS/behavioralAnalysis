function chance = empiricalChancePerf(lg)

% chance = empiricalChancePerf(lg)
% calculates empirical chance performance based on actual trial draws

mice    = unique(lg.mouseID);
nmice   = numel(mice);
ntrials = nan(nmice,1);
frac    = ntrials;

% empirical trial fraction per mouse
for iMouse = 1:nmice
  trialType       = lg.trialType(lg.mouseID == mice(iMouse));
  frac(iMouse)    = sum(trialType == 1)./numel(trialType);
  ntrials(iMouse) = numel(trialType);
end

% weighted sume across mice
frac(frac>.5) = 1-frac(frac>.5);
chance        = sum(frac.*ntrials)/sum(ntrials); 