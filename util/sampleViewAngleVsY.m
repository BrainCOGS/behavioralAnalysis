% [viewAngle, yPos, yIndex] = sampleViewAngleVsY(position, ySample, doFirst)
% Computes the view angle as a function of y-position in the maze, evaluated at the given ySample
% values. Since this is not necessarily unique (e.g. the mouse can double back), first we compute
% the cumulative maximum y and then define view angle as sampled using those y values. Specifically
% the last value of y that is greater than or equal to ySample is selected
% (unless a 3rd input logical true is passed, in which case it does first)
%
% Except for yPos which is a cell array of different lengths (depending on the duration of the
% trial), this function returns n-by-t matrices where n is the number of ySample values and t is the
% number of trials. 
% 
% by Sue Ann Koay

function [viewAngle, yPos, yIndex] = sampleViewAngleVsY(position, ySample, doFirst)

  if nargin < 3; doFirst = false; end % sample 1st instead of last unique y
  if doFirst; warning('deprecated code, only last unique position is supported now'); end
  
  yPos        = cellfun(@(x) cummax(x(:,2)), position, 'UniformOutput', false);
  yIndex      = accumfun(2, @(x) binarySearch(x, ySample, 1, 1)', yPos);
 
  viewAngle = accumfun(2, @(x) position{x}(yIndex(:,x),end), 1:numel(position));
  
end
