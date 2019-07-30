function lsrPerf = mouseLsrStats(lsrPerf,alpha)

% lsrPerf = mouseLsrStats(lsrPerf,alpha)
% signed rank tests (or paired t where adequate) across mice for laser
% for each location will create structure mouseStats with descriptive stats
% and pvalues for
% lsrPerf is output structure of analyzeConcatLsrLog()
% alpha is significance level
% indicators = {'percCorrect','percCorrect_R','percCorrect_L','bias','speed',...
%              'percFinish','excessTravel','dur','viewOverlap'};

if nargin < 2; alpha = 0.05; end
if ~iscell(lsrPerf); lsrPerf = {lsrPerf}; revertBack = true; else revertBack = false; end

indicators = {'percCorrect','percCorrect_R','percCorrect_L','percCorrect_easy'  ...
             ,'percCorrect_hard','bias','bias_abs','lapse','percFinish'         ...
             ,'excessTravel','dur','speed','unusualEvents'};

for iVec = 1:numel(indicators)
  if ~isfield(lsrPerf{1}.ctrl,[indicators{iVec} '_mouse']); continue; end
  for iLoc = 1:numel(lsrPerf)
    thisctrl     = lsrPerf{iLoc}.ctrl.([indicators{iVec} '_mouse']);
    thislsr      = lsrPerf{iLoc}.lsr.([indicators{iVec} '_mouse']);
    [p,testname] = calculatePValue(thislsr,thisctrl,'classicPaired',.05);
    lsrPerf{iLoc}.stats.xMouse.(['p_' indicators{iVec}])        = double(p);
    lsrPerf{iLoc}.stats.xMouse.(['testname_' indicators{iVec}]) = testname;
    
    if isfield(lsrPerf{iLoc},'cueHalf')
      thisctrl     = lsrPerf{iLoc}.cueHalf.ctrl.([indicators{iVec} '_mouse']);
      thislsr      = lsrPerf{iLoc}.cueHalf.lsr.([indicators{iVec} '_mouse']);
      [p,testname] = calculatePValue(thislsr,thisctrl,'classicPaired',.05);
      lsrPerf{iLoc}.cueHalf.stats.xMouse.(['p_' indicators{iVec}])        = double(p);
      lsrPerf{iLoc}.cueHalf.stats.xMouse.(['testname_' indicators{iVec}]) = testname;
    end
    if isfield(lsrPerf{iLoc},'cueHalf_vsCtrl')
      thisctrl     = lsrPerf{iLoc}.cueHalf_vsCtrl.ctrl.([indicators{iVec} '_mouse']);
      thislsr      = lsrPerf{iLoc}.cueHalf_vsCtrl.lsr.([indicators{iVec} '_mouse']);
      [p,testname] = calculatePValue(thislsr,thisctrl,'classicPaired',.05);
      lsrPerf{iLoc}.cueHalf_vsCtrl.stats.xMouse.(['p_' indicators{iVec}])        = double(p);
      lsrPerf{iLoc}.cueHalf_vsCtrl.stats.xMouse.(['testname_' indicators{iVec}]) = testname;
    end
  end
  % FDR correction
  if iLoc > 1 
    pvals                  = cellfun(@(x)(x.stats.xMouse.(['p_' indicators{iVec}])),lsrPerf);
    [lsrPerf{1}.stats.xMouse.(['isSig_' indicators{iVec}]), lsrPerf{1}.stats.xMouse.(['alpha_correct_' indicators{iVec}])] ...
                           = FDR(pvals', alpha);
    if isfield(lsrPerf{1},'cueHalf')
      pvals                  = cellfun(@(x)(x.cueHalf.stats.xMouse.(['p_' indicators{iVec}])),lsrPerf);
      [lsrPerf{1}.cueHalf.stats.xMouse.(['isSig_' indicators{iVec}]), lsrPerf{1}.cueHalf.stats.xMouse.(['alpha_correct_' indicators{iVec}])] ...
                           = FDR(pvals', alpha);
    end
    if isfield(lsrPerf{1},'cueHalf_vsCtrl')
      pvals                  = cellfun(@(x)(x.cueHalf_vsCtrl.stats.xMouse.(['p_' indicators{iVec}])),lsrPerf);
      [lsrPerf{1}.cueHalf_vsCtrl.stats.xMouse.(['isSig_' indicators{iVec}]), lsrPerf{1}.cueHalf_vsCtrl.stats.xMouse.(['alpha_correct_' indicators{iVec}])] ...
                           = FDR(pvals', alpha);
    end
  end
end

if revertBack; lsrPerf = lsrPerf{1}; end
