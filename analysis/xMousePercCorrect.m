function pc = xMousePercCorrect(lg,minNumTrials,doPercGoodSess)

% pc = xMousePercCorrect(lg,minNumTrials,doPercGoodSess)
% average per-mouse % correct performance
% lg is flattened behavior log
% minNumTrials is minimal trial count to be included in analysis (default
% 1000)
% doPercGoodSess is boolean flag to calculate fracton of sessions with good
% blocks (> 60% correct)

if nargin < 2; minNumTrials   = 1000; end
if nargin < 3; doPercGoodSess = true; end

lg    = cleanupConcatLog(lg, minNumTrials);

% overall with binomial CIs
stdInterval      = normcdf(1, 0, 1) - normcdf(-1, 0, 1);
pc.overall       = sum(lg.trialType == lg.choice)./length(lg.trialType)*100;
[~,pc.overallCI] = binofit(sum(lg.trialType == lg.choice),length(lg.trialType),1-stdInterval);

% by mouse
pc.mice         = unique(lg.mouseID);
pc.nmice        = length(pc.mice);
pc.mousePC      = zeros(1,pc.nmice);
pc.goodSessPerc = nan(1,pc.nmice);

for iMouse = 1:pc.nmice
  [trialType,choice] = selectMouseTrials(lg, pc.mice(iMouse), 'trialType', 'choice');
  pc.mousePC(iMouse) = sum(trialType == choice)/numel(choice)*100;
  
  % get file list to estimate number of good sessions
  if doPercGoodSess
    filters.mazetype        = 'last';
    filters.mouseID         = analysisParams.miceBehav(pc.mice(iMouse));
    fl                      = filterLogFilesBehav(filters);
    total                   = numel(fl);
    good                    = numel(unique(lg.sessionID(lg.mouseID == pc.mice(iMouse))));
    pc.goodSessPerc(iMouse) = 100*(good/total);
  end
end

pc.mouseAvg         = mean(pc.mousePC);
pc.mouseSem         = std(pc.mousePC)./sqrt(pc.nmice-1);
pc.mouseAvgGoodSess = mean(pc.goodSessPerc);
pc.mouseSemGoodSess = std(pc.goodSessPerc)./sqrt(pc.nmice-1);