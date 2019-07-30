function [p,testname,sd] = calculatePValue(vec,vecCtrl,method,alpha)

% [p,testname,sd] = calculatePValue(vec,vecCtrl,method,alpha)
%   calculates significance of difference between two vectors, vc and
%   vecCtrl using the method as defined below
%   method: 'boot' for bootstrapping, 'classic' for ranksum or ttest
%   according to normality test using lillietest, 'classicPaired' for
%   paired equivalents of classic
%   alpha is significance threshold
%   
% p is pvalue
% testname is tring with test used
% sd is standrad deviations of bootstrapping iterations if applicable

%%  initialize
if nargin < 2
  vecCtrl = zeros(size(vec));
end
if nargin < 3
  method  = 'boot';
end
if nargin < 4
  alpha   = 0.01;
end

if isempty(vecCtrl); vecCtrl = zeros(size(vec)); end

if size(vec,1) < size(vec,2)
  vec     = vec';
  vecCtrl = vecCtrl';
end

switch method
  case 'boot'
    % generic bootstrapping (two-tailed)
    testname = 'boot';
    niter    = 10000;
    bootvec  = zeros(1,niter);
    vecAll   = [vec; vecCtrl];
    vecTags  = [true(size(vec)); false(size(vecCtrl))];
    for ii = 1:niter
      idx         = randsample(length(vecAll),length(vecAll),true); % sample with replacement
      ivecAll     = vecAll(idx);
      ivecTags    = vecTags(idx);
      bootvec(ii) = mean(ivecAll(ivecTags)) - mean(ivecAll(~ivecTags));
    end
    p         = sum(sign(bootvec) ~= sign(mean(vec) - mean(vecCtrl))) / niter;
    sd        = std(bootvec);
    
  case 'classic'
    % test for normality and choose appropriate method
    sd                  = [];
    nanEntries          = isnan(vec);
    vec(nanEntries)     = [];
    nanEntries          = isnan(vecCtrl);
    vecCtrl(nanEntries) = [];
    try
    if lillietest(vec) || (sum(vecCtrl>0 & lillietest(vecCtrl)))
      p        = ranksum(vec,vecCtrl,alpha);
      testname = 'ranksum';
    else
      [~,p]    = ttest2(vec,vecCtrl,alpha);
      testname = 'ttest';
    end
    catch
      p = nan;
      testname = '';
    end
  case 'classicPaired'
    % test for normality and choose appropriate method
    sd                  = [];
    nanEntries          = isnan(vec) | isnan(vecCtrl);
    vec(nanEntries)     = [];
    vecCtrl(nanEntries) = [];
    try
    try
      if lillietest(vec) || (sum(vecCtrl>0 & lillietest(vecCtrl)))
        p        = signrank(vec,vecCtrl,alpha);
        testname = 'signrank';
      else
        [~,p]    = ttest(vec,vecCtrl,alpha);
        testname = 'ttest';
      end
    catch
      p        = signrank(vec,vecCtrl,alpha);
      testname = 'signrank';
    end
    catch
      p = nan;
      testname = '';
    end
end