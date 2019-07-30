% Finds trials with at least that nHistory trials past in the log that consecutively lead up to
% trials in the returned present vector.
%
% Input:
%   nPast_good  : from cleanupConcatLog()
% Output:
%   present     : vector of indices of trials with at least nHistory past trials
%   past        : for convenience, each column i of this is the index of the i-th back history
%                 trial, i.e. present - i
function [present, past] = findTrialsWithHistory(nPast_good, nHistory)
  
  present   = find(nPast_good(:) >= nHistory);
  if nargout > 1
    past    = bsxfun(@minus, present, 1:nHistory);
  end
  
end
