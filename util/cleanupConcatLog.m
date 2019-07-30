function [newlog, info] = cleanupConcatLog(lg, minNT, metaMouse, minHistory, selectTrialsWithHistory, relaxCriteria, deleteCueHalfTrials, perfTh)

% [newlog, info] = cleanupConcatLog(lg, minNT, metaMouse, minHistory, selectTrialsWithHistory, relaxCriteria, deleteCueHalfTrials, perfTh)
% apply trial selections in isBadTrial and select mice with at least minNT
% trials. If isempty(minNT) all mice will be selected. Output is new log
% structure containing only relevant trials
%
% minHistory > 0 will apply the minNT selection on the number of trials with
% that many consecutive trials past, as defined by the findTrialsWithHistory() 
% function. A field 'nPast_good' is added to newlog to be used for this selection.
%
% selectTrialsWithHistory = true should be used to remove trials without at
% least minHistory trials consecutively in the past, from the resulting newlog.
% By default it is false; note that trials without a sufficient number of
% history trials should be retained in the logs if one wants to select n-back
% trials at some point in the analysis.
%
% relaxCriteria = true will relax the selection criteria in isBadTrial, in
% particular the excessTravel requirement is removed. This should be used
% for trial history studies in order to avoid biases due to strict selections.
% If not provided, it will be assumed to be false.
%
% deleteCueHalfTrials = true will remove special control trials present in optogenetic
% inactivation dataset where all cues are deleted from either the first or last 100 cm.
% If not provided, it will be assumed to be false.

if nargin < 2
   minNT      = 1000;
end
if nargin < 3
  metaMouse   = false;
end
if nargin < 4 || isempty(minHistory)
  minHistory = 0;
end
if nargin < 5 || isempty(selectTrialsWithHistory)
  selectTrialsWithHistory = false;
end
if nargin < 6 || isempty(relaxCriteria)
%   relaxCriteria = minHistory > 0;
  relaxCriteria = false;
end
if nargin < 7 || isempty(deleteCueHalfTrials)
  deleteCueHalfTrials = false;
end
if nargin < 8 || isempty(perfTh)
  perfTh = .6;
end
%% add to the log a counter for the number of valid historical trials
badtrial = isBadTrial(lg, [], relaxCriteria, perfTh);
nHistory = nan(size(lg.trialType));
goodCount = 0;
for iTrial = 1:numel(nHistory)
  % First trial of a block has no history
  if lg.firstTrialofBlock(iTrial)
    goodCount = 0;
  end
  nHistory(iTrial)  = goodCount;
  
  % Bad trials have a history but break the consecutive counts
  if badtrial(iTrial)
    goodCount = 0;
  else
    goodCount = goodCount + 1;
  end
end
lg.nPast_good = nHistory;

%% if desired exclude special laser control trials with no cues in either 2st or 2nd halves
if deleteCueHalfTrials && isfield(lg,'laserON')
  isCtrlTrial = (cellfun(@(x,y)(sum(x > 100) == 0 & sum(y > 100) == 0),lg.cuePos_R,lg.cuePos_L) | ...
                cellfun(@(x,y)(sum(x < 100) == 0 & sum(y < 100) == 0),lg.cuePos_R,lg.cuePos_L)) & ...
                (lg.nCues_RminusL == 0 | (lg.nCues_RminusL > 0 & lg.trialType == analysisParams.leftCode) | ...
                (lg.nCues_RminusL < 0 & lg.trialType == analysisParams.rightCode)) & ~lg.laserON;
  badtrial    = badtrial | isCtrlTrial;
end

%% veto trials without a sufficient number of consecutive history trials, if so desired
if selectTrialsWithHistory
  badtrial  = badtrial | lg.nPast_good < minHistory;
end


%% select out bad trials first
fls      = fieldnames(lg);
newlog   = struct();
for iF = 1:numel(fls)
   if ~strcmpi(fls{iF},'keyFrameLabels') && ~strcmpi(fls{iF},'mouseNames')
      newlog.(fls{iF}) = lg.(fls{iF})(~badtrial);
   end
end


%% apply number of trials criterion
mice       = unique(newlog.mouseID);
lowtrialct = false(size(newlog.mouseID));
if ~isempty(minNT)
   for iMouse = 1:numel(mice)
      sel = (newlog.mouseID == mice(iMouse));
      if sum(sel & newlog.nPast_good >= minHistory) < minNT
        lowtrialct(sel) = true; 
      end
   end
end

for iF = 1:numel(fls)
   if ~strcmpi(fls{iF},'keyFrameLabels') && ~strcmpi(fls{iF},'mouseNames')
      newlog.(fls{iF}) = newlog.(fls{iF})(~lowtrialct);
   end
end

%% update info
mice                  = unique(newlog.mouseID);
newlog.mouseNames     = analysisParams.miceBehav(mice);

if isfield(lg,'keyFrameLabels')
  newlog.keyFrameLabels = lg.keyFrameLabels;
end

if metaMouse
  mice(end+1)     = -1;
  newlog.mouseNames{end+1}  = 'metamouse';
end

info              = struct();
info.minNumTrials = minNT;
info.mouseIndices = mice;
info.mouseNames   = newlog.mouseNames;

%% Addtional selection info for convenience and bookkeeping
info.mouseNTrials = arrayfun(@(x) sum(newlog.mouseID == x & newlog.nPast_good >= minHistory), mice);
info.totalTrials  = sum(info.mouseNTrials);
if metaMouse
  info.mouseNTrials(end)  = info.totalTrials;
end
