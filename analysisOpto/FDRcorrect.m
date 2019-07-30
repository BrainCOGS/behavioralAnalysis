function [alpha_correct,isSig] = FDRcorrect(lsrPerf,alpha,whichP)

% [alpha_correct,isSig] = FDRcorrect(lsrPerf,alpha,whichP)
%
% performs a false discovery rate correction according to method described
% in Benjamini & Hochberg, 1995, J Royal Stat Soci B, 57(1):289-300 and
% used in Guo et al (2014) Neuron
% Briefly, it ranks p-values in ascending order and defines a value i as being
% significant if it satisfies P(i) <= alpha*i/n where n is the number of
% comparisons
% lsrPerf is cell array containing data structures with stats for each
% location, alpha is orginal value (def. 0.025 assuming two-tailed tests),
% whichP is string to indicate which behavioral parameter (see function)
% alpha_correct is FDR-corrected statistical threshold
% isSig is nlocations x 1 boolean vector where true indicates significance
%
% LP aug 2016

if nargin < 2
  alpha = analysisParams.bootalpha;
end
if nargin < 3
  whichP = 'percCorrect';
end

n     = length(lsrPerf);
pvals = zeros(1,n);
switch whichP
  case 'percCorrect'
    for ii = 1:n; pvals(ii) = lsrPerf{ii}.stats.bootPerf.p_percCorrect; end
  case 'bias'
    for ii = 1:n; pvals(ii) = lsrPerf{ii}.stats.bootPerf.p_bias; end
  case 'bias_abs'
    for ii = 1:n; pvals(ii) = lsrPerf{ii}.stats.bootPerf.p_bias_abs; end
  case 'lapse'
    for ii = 1:n; pvals(ii) = lsrPerf{ii}.stats.bootPerf.p_lapse; end
  case 'percCorrect_easy'
    for ii = 1:n; pvals(ii) = lsrPerf{ii}.stats.bootPerf.p_percCorrect_easy; end
  case 'percCorrect_hard'
    for ii = 1:n; pvals(ii) = lsrPerf{ii}.stats.bootPerf.p_percCorrect_hard; end
  case 'percCorrect_R'
    for ii = 1:n; pvals(ii) = lsrPerf{ii}.stats.bootPerf.p_percCorrect_R; end
  case 'percCorrect_L'
    for ii = 1:n; pvals(ii) = lsrPerf{ii}.stats.bootPerf.p_percCorrect_L; end
  case 'percFinish'
    for ii = 1:n; pvals(ii) = lsrPerf{ii}.stats.bootPerf.p_percFinish; end
  case 'speed'
    for ii = 1:n; pvals(ii) = lsrPerf{ii}.stats.p_speed; end
  case 'dur'
    for ii = 1:n; pvals(ii) = lsrPerf{ii}.stats.p_dur; end
  case 'excessTravel'
    for ii = 1:n; pvals(ii) = lsrPerf{ii}.stats.p_excessTravel; end
  case 'unusualEvents'
    for ii = 1:n; pvals(ii) = lsrPerf{ii}.stats.bootMotor.p_unusualEvents; end
  case 'slope'
    for ii = 1:n; pvals(ii) = lsrPerf{ii}.stats.bootPerf.p_slope; end
  case 'mousePcorrect'
    for ii = 1:n; pvals(ii) = lsrPerf{ii}.mouseStats.lsr.p_percCorrect; end
  case 'viewDecode'
    for ii = 1:n; pvals(ii) = lsrPerf{ii}.viewAngle.p_decodeAcc_overall; end
  case 'viewAngSD'
    for ii = 1:n; pvals(ii) = lsrPerf{ii}.viewAngle.p_viewAngSD; end
  case 'logRegDecayRatio'
    for ii = 1:n; pvals(ii) = lsrPerf{ii}.lsr.logisticReg.decayIndex_p; end
  case 'logRegAvgWeight'
    for ii = 1:n; pvals(ii) = lsrPerf{ii}.lsr.logisticReg.avgWeight_p; end
%   case 'viewOverlap'
%     for ii = 1:n; pvals(ii) = lsrPerf{ii}.stats.bootPerf.p_viewOverlap; end
end

[isSig, alpha_correct] = FDR(pvals', alpha);
isSig                  = isSig';