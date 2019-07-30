function [nr,nc] = subplotOrg(nplots,maxcols)

%[nr,nc] = subplotOrg(nplots,maxcols)
%returns numbers of rows and columns in a multi-panel figure; given the
%number of panels and the maximum desired number of columns

if nargin < 2
  maxcols = 3;
end

nr = ceil(nplots/maxcols);
nc = ceil(nplots/nr);