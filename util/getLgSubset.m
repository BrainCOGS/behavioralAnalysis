function sublg = getLgSubset(lg,trialIdx)

% sublg = getLgSubset(lg,trialIdx)
% generate sublg with subset of trials trialIdx from lg


varls   = fieldnames(lg);
nTrials = numel(lg.choice);
for iVar = 1:numel(varls)
  if numel(lg.(varls{iVar})) < nTrials
    sublg.(varls{iVar}) = lg.(varls{iVar});
  else
    sublg.(varls{iVar}) = lg.(varls{iVar})(trialIdx);
  end
end
