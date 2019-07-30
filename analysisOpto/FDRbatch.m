function lsrPerf = FDRbatch(lsrPerf,alpha,indicators)

% lsrPerf = FDRbatch(lsrPerf,alpha,indicators)
% performs FDR correction on a list of indicators
% lsrPerf is output structure of analyzeConcatLsrLog()
% alpha is significance level before correction
% indicator is cell array with list of indicators (strings)

if nargin < 2 || isempty(alpha)
  alpha = analysisParams.bootalpha;
end
if nargin < 3
  indicators = {'percCorrect','percCorrect_L','percCorrect_R','bias',               ...
    'speed','dur','percFinish','excessTravel','slope','viewDecode','unusualEvents', ...
    'percCorrect_easy','percCorrect_hard','bias_abs','lapse','viewAngSD',           ...
    'logRegDecayRatio','logRegAvgWeight'};
end
for ii = 1:length(indicators)
%   try
    ac = FDRcorrect(lsrPerf,alpha,indicators{ii});
%   catch
%     ac = 0.025;
%   end
  eval(sprintf('lsrPerf{1}.FDR.alpha_correct_%s = ac;',indicators{ii}))
end

lsrPerf{1}.alpha = alpha;