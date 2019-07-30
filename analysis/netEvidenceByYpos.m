function [nev,yPerEv,bins,evs] = netEvidenceByYpos(lg)

% [nev,yPerEv,bins,evs] = netEvidenceByYpos(lg)
% extracts net evidence per y position for each trial in flattened log, lg

bins    = 10:10:200;
evs     = 0:15;
choice  = lg.choice(~isnan(lg.choice) & lg.choice ~= -1);
et      = lg.excessTravel(~isnan(lg.choice) & lg.choice ~= -1);
cr      = lg.cuePos_R(~isnan(lg.choice) & lg.choice ~= -1);
cl      = lg.cuePos_L(~isnan(lg.choice) & lg.choice ~= -1);
ntrials = numel(choice);
nev     = zeros(ntrials,numel(bins));
yPerEv  = nan(ntrials,numel(evs));

for ii = 1:ntrials
  if et(ii) > .1; continue; end
  for jj = 1:numel(bins)
    nR         = sum(cr{ii} <= bins(jj));
    nL         = sum(cl{ii} <= bins(jj));
    nev(ii,jj) = abs(nR-nL);
  end
  for jj = 1:numel(evs)
    temp = find(nev(ii,:) >= evs(jj),1,'first');
    if ~isempty(temp)
      yPerEv(ii,jj) = bins(temp);
    end
  end
end