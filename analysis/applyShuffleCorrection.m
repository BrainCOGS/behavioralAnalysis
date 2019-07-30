function lsrPerf = applyShuffleCorrection(lsrPerf,lsrPerfShuffle,alpha)

% lsrPerf = applyShuffleCorrection(lsrPerf,lsrPerfShuffle,alpha)
% pick p value from shuffled data such that alpha% of regions have
% significant effects (ie do heuristic multiple comparisons correction)

if nargin < 3; alpha = analysisParams.bootalpha; end

%%
indicators = {'percCorrect','percCorrect_L','percCorrect_R','bias',...
  'speed','dur','percFinish','excessTravel','slope','viewDecode','unusualEvents',...
  'percCorrect_easy','percCorrect_hard','bias_abs','lapse','logRegSlope','logRegDecayIndex'};
nReg       = round(numel(lsrPerf)*alpha);

%%
for iVar = 1:length(indicators)
  pvals = min(getPvals(lsrPerfShuffle,indicators{iVar})); % min pval for each region across shuffles
  pvals = rank(pvals,'ascend');
  ac    = pvals(nReg+1);
  lsrPerf{1}.FDR.(['alpha_correct_' indicators{iVar}]) = ac;
end

lsrPerf{1}.alpha = alpha;

end

%% get appropriate pval
function pvals = getPvals(lsrPerf,whichP)

nLoc  = numel(lsrPerf);
nShuf = numel(lsrPerf{1});
pvals = zeros(nShuf,nLoc);
switch whichP
  case 'percCorrect'
    for iSh = 1:nShuf; for iLoc = 1:nLoc; pvals(iSh,iLoc) = lsrPerf{iLoc}{iSh}.stats.bootPerf.p_percCorrect; end; end
  case 'bias'
    for iSh = 1:nShuf; for iLoc = 1:nLoc; pvals(iSh,iLoc) = lsrPerf{iLoc}{iSh}.stats.bootPerf.p_bias; end; end
  case 'bias_abs'
    for iSh = 1:nShuf; for iLoc = 1:nLoc; pvals(iSh,iLoc) = lsrPerf{iLoc}{iSh}.stats.bootPerf.p_bias_abs; end; end
  case 'lapse'
    for iSh = 1:nShuf; for iLoc = 1:nLoc; pvals(iSh,iLoc) = lsrPerf{iLoc}{iSh}.stats.bootPerf.p_lapse; end; end
  case 'percCorrect_easy'
    for iSh = 1:nShuf; for iLoc = 1:nLoc; pvals(iSh,iLoc) = lsrPerf{iLoc}{iSh}.stats.bootPerf.p_percCorrect_easy; end; end
  case 'percCorrect_hard'
    for iSh = 1:nShuf; for iLoc = 1:nLoc; pvals(iSh,iLoc) = lsrPerf{iLoc}{iSh}.stats.bootPerf.p_percCorrect_hard; end; end
  case 'percCorrect_R'
    for iSh = 1:nShuf; for iLoc = 1:nLoc; pvals(iSh,iLoc) = lsrPerf{iLoc}{iSh}.stats.bootPerf.p_percCorrect_R; end; end
  case 'percCorrect_L'
    for iSh = 1:nShuf; for iLoc = 1:nLoc; pvals(iSh,iLoc) = lsrPerf{iLoc}{iSh}.stats.bootPerf.p_percCorrect_L; end; end
  case 'percFinish'
    for iSh = 1:nShuf; for iLoc = 1:nLoc; pvals(iSh,iLoc) = lsrPerf{iLoc}{iSh}.stats.bootPerf.p_percFinish; end; end
  case 'speed'
    for iSh = 1:nShuf; for iLoc = 1:nLoc; pvals(iSh,iLoc) = lsrPerf{iLoc}{iSh}.stats.p_speed; end; end
  case 'dur'
    for iSh = 1:nShuf; for iLoc = 1:nLoc; pvals(iSh,iLoc) = lsrPerf{iLoc}{iSh}.stats.p_dur; end; end
  case 'excessTravel'
    for iSh = 1:nShuf; for iLoc = 1:nLoc; pvals(iSh,iLoc) = lsrPerf{iLoc}{iSh}.stats.p_excessTravel; end; end
  case 'unusualEvents'
    for iSh = 1:nShuf; for iLoc = 1:nLoc; pvals(iSh,iLoc) = lsrPerf{iLoc}{iSh}.stats.bootMotor.p_unusualEvents; end; end
  case 'slope'
    for iSh = 1:nShuf; for iLoc = 1:nLoc; pvals(iSh,iLoc) = lsrPerf{iLoc}{iSh}.stats.bootPerf.p_slope; end; end
  case 'mousePcorrect'
    for iSh = 1:nShuf; for iLoc = 1:nLoc; pvals(iSh,iLoc) = lsrPerf{iLoc}{iSh}.mouseStats.lsr.p_percCorrect; end; end
  case 'viewDecode'
    for iSh = 1:nShuf; for iLoc = 1:nLoc; pvals(iSh,iLoc) = lsrPerf{iLoc}{iSh}.viewAngle.p_decodeAcc_overall; end; end
  case 'logRegSlope'
    for iSh = 1:nShuf; for iLoc = 1:nLoc; pvals(iSh,iLoc) = lsrPerf{iLoc}{iSh}.lsr.logisticReg.slope_p; end; end
  case 'logRegDecayIndex'
    for iSh = 1:nShuf; for iLoc = 1:nLoc; pvals(iSh,iLoc) = lsrPerf{iLoc}{iSh}.lsr.logisticReg.decayIndex_p; end; end
%   case 'viewOverlap'
%     for iSh = 1:nShuf; for iLoc = 1:nLoc; pvals(iSh,iLoc) = lsrPerf{iLoc}{iSh}.stats.bootPerf.p_viewOverlap; end; end
end

end