function analyzeConcatLsrLog_batch(exptType,updateDataset,diagnoseDataset,minNumTrialsPerMouse)

% analyzeConcatLsrLog_batch(exptType,updateDataset,diagnoseDataset,minNumTrialsPerMouse)
% compiles relevant experiments and analyzes inactivation effects
% exptType: 'wholeTrial6mw', 'wholeTrial1mw','wholeTrial2mw','wholeTrialClust', 
%           'WTctrls', 'epochs', 'memGuide', 'visGuide'

%% initialize
if nargin < 1; exptType        = 'wholeTrial6mw'; end
if nargin < 2; updateDataset        = false;  end
if nargin < 3; diagnoseDataset      = false;  end
if nargin < 4; minNumTrialsPerMouse = 0;      end

spockFlag = isThisSpock;
poolobj = gcp('nocreate');
if isempty(poolobj); poolobj = parpool; end

ap                       = analysisParams;
cfg.nBins_logReg         = 4;
cfg.minNumTrialsPerMouse = minNumTrialsPerMouse;

try
%% update logs, save data diagnosis
if updateDataset;   summarizeLsrExpt_batch; end
if diagnoseDataset; diagnoseLsrDataset;     end
warning('off','all'); 

switch exptType
  % full grid whole trial, 6 mW, vgat
  case 'wholeTrial6mw' 
    %% full grid whole trial, 6 mW, vgat

    thist = tic;
    fprintf('-------------------\n\n')

    % get grid locations, labels
    filters.mazetype                = {'last';'last-1'};
    filters.perfTh                  = 0.6;
    filters.deleteNonMainMazeTrials = true;
    filters.powerLogic              = '==';
    filters.power                   = 6;
    filters.epoch                   = 'whole';
    filters.genotype                = 'vgat';
    filters.gridLabel               = {'fullGridBilateral.mat';      ...
                                       'fullLessGridBilateral.mat';  ...     
                                       'posteriorGridBilateral.mat'; ...
                                       'anteriorExtraGridBilateral.mat'};

    lg               = concatLogsLsr([],filters,spockFlag);
    if iscell(filters.mazetype); mazetype = filters.mazetype{1}; else mazetype = filters.mazetype; end
    fn               = sprintf('fullGrid_%sMaze_%s_%s_perfTh%d',mazetype,filters.epoch,filters.genotype,100*filters.perfTh);
    if isfield(filters,'power'); fn = sprintf('%s_%dmW',fn,round(filters.power)); end
    fn               = sprintf('%s_nBinLogReg%d',fn,cfg.nBins_logReg);

    dataLbl          = sprintf('fullGrid, maze: %s, epoch: %s, %s mice, perfTh: %d%%',mazetype,filters.epoch,filters.genotype,100*filters.perfTh);
    if isfield(filters,'power'); dataLbl = sprintf('%s, %dmW',dataLbl,round(filters.power)); end

    data             = load('fullGridBilateral.mat','grid');
    grid             = data.grid;
    locIdx           = getMasterLocationIdx(grid);

    % loop over locations
    parfor iLoc = 1:numel(locIdx)
      lsrPerf{iLoc} = analyzeConcatLsrLog(lg,locIdx(iLoc),filters,cfg,grid);
    end
    lsrPerf         = FDRbatch(lsrPerf);
    lsrPerf         = mouseLsrStats(lsrPerf);

    % plot and save
    rootdir         = ap.getRootDir(spockFlag);
    cd(rootdir)
    save(fn,'lsrPerf','filters','fn','dataLbl')

    fprintf('\nANALYZED %s in %1.1f min\n',upper(fn),toc(thist)/60)
    clear grid lsrPerf locIdx filters

    %% full grid whole trial, 2 mW, vgat
    case 'wholeTrial2mw' 
      
      thist = tic;
      fprintf('-------------------\n\n')

      filters.mazetype                = {'last';'last-1'};
      filters.perfTh                  = 0.6;
      filters.deleteNonMainMazeTrials = true;
      filters.epoch                   = 'whole';
      filters.genotype                = 'vgat';
      filters.gridLabel               = {'fullGridBilateral.mat';      ...
                                         'fullLessGridBilateral.mat';  ...     
                                         'posteriorGridBilateral.mat'; ...
                                         'anteriorExtraGridBilateral.mat'};
      filters.powerLogic              = '==';
      filters.power                   = 2;

      % get grid locations, labels
      lg               = concatLogsLsr([],filters,spockFlag);
      if iscell(filters.mazetype); mazetype = filters.mazetype{1}; else mazetype = filters.mazetype; end
      fn               = sprintf('fullGrid_%sMaze_%s_%s_perfTh%d',mazetype,filters.epoch,filters.genotype,100*filters.perfTh);
      if isfield(filters,'power'); fn = sprintf('%s_%dmW',fn,round(filters.power)); end
      fn               = sprintf('%s_nBinLogReg%d',fn,cfg.nBins_logReg);

      dataLbl          = sprintf('fullGrid, maze: %s, epoch: %s, %s mice, perfTh: %d%%',mazetype,filters.epoch,filters.genotype,100*filters.perfTh);
      if isfield(filters,'power'); dataLbl = sprintf('%s, %dmW',dataLbl,round(filters.power)); end

      data             = load('fullGridBilateral.mat','grid');
      grid             = data.grid;
      locIdx           = getMasterLocationIdx(grid);

      % loop over locations
      parfor iLoc = 1:numel(locIdx)
        lsrPerf{iLoc} = analyzeConcatLsrLog(lg,locIdx(iLoc),filters,cfg,grid);
      end
      lsrPerf         = FDRbatch(lsrPerf);
      lsrPerf         = mouseLsrStats(lsrPerf);

      % plot and save
      rootdir         = ap.getRootDir(spockFlag);
      cd(rootdir)
      save(fn,'lsrPerf','filters','fn','dataLbl')

      fprintf('\nANALYZED %s in %1.1f min\n',upper(fn),toc(thist)/60)
      clear filters grid lsrPerf locIdx
      
    %% full grid whole trial, 1 mW, vgat
    case 'wholeTrial1mw' 
      
      thist = tic;
      fprintf('-------------------\n\n')

      filters.mazetype                = {'last';'last-1'};
      filters.perfTh                  = 0.6;
      filters.deleteNonMainMazeTrials = true;
      filters.epoch                   = 'whole';
      filters.genotype                = 'vgat';
      filters.gridLabel               = {'fullGridBilateral.mat';      ...
                                         'fullLessGridBilateral.mat';  ...     
                                         'posteriorGridBilateral.mat'; ...
                                         'anteriorExtraGridBilateral.mat'};
      filters.powerLogic              = '==';
      filters.power                   = 1;

      % get grid locations, labels
      lg               = concatLogsLsr([],filters,spockFlag);
      if iscell(filters.mazetype); mazetype = filters.mazetype{1}; else mazetype = filters.mazetype; end
      fn               = sprintf('fullGrid_%sMaze_%s_%s_perfTh%d',mazetype,filters.epoch,filters.genotype,100*filters.perfTh);
      if isfield(filters,'power'); fn = sprintf('%s_%dmW',fn,round(filters.power)); end
      fn               = sprintf('%s_nBinLogReg%d',fn,cfg.nBins_logReg);
%       fn               = sprintf('%s_minNumTrials%d',fn,cfg.minNumTrialsPerMouse);
      
      dataLbl          = sprintf('fullGrid, maze: %s, epoch: %s, %s mice, perfTh: %d%%',mazetype,filters.epoch,filters.genotype,100*filters.perfTh);
      if isfield(filters,'power'); dataLbl = sprintf('%s, %dmW',dataLbl,round(filters.power)); end

      data             = load('fullGridBilateral.mat','grid');
      grid             = data.grid;
      locIdx           = getMasterLocationIdx(grid);

      % loop over locations
      parfor iLoc = 1:numel(locIdx)
        lsrPerf{iLoc} = analyzeConcatLsrLog(lg,locIdx(iLoc),filters,cfg,grid);
      end
      lsrPerf         = FDRbatch(lsrPerf);
      lsrPerf         = mouseLsrStats(lsrPerf);

      % plot and save
      rootdir         = ap.getRootDir(spockFlag);
      cd(rootdir)
      save(fn,'lsrPerf','filters','fn','dataLbl')

      fprintf('\nANALYZED %s in %1.1f min\n',upper(fn),toc(thist)/60)
      clear filters grid lsrPerf locIdx

  %% combine locations from full grid according to clustering
  case 'wholeTrialClust'
    thist = tic;
    fprintf('-------------------\n\n')

    % get grid locations, labels
    filters.mazetype                = {'last';'last-1'};
    filters.perfTh                  = 0.6;
    filters.deleteNonMainMazeTrials = true;
    filters.powerLogic              = '==';
    filters.power                   = 6;
    filters.epoch                   = 'whole';
    filters.genotype                = 'vgat';
    filters.gridLabel               = {'fullGridBilateral.mat';      ...
                                       'fullLessGridBilateral.mat';  ...     
                                       'posteriorGridBilateral.mat'; ...
                                       'anteriorExtraGridBilateral.mat'};

    lg        = concatLogsLsr([],filters,spockFlag);            

    loclbls   = {'clust1','clust2','clust3'};
    coords{1} = {[-2 -2.5; 2 -2.5];[-2 -3.5; 2 -3.5];[-3 -3.5; 3 -3.5];[-3 -2.5; 3 -2.5]; ...
                 [-2 -1.5; 2 -1.5];[-3 -1.5; 3 -1.5];[-1 -3.5; 1 -3.5];[-1 -2.5; 1 -2.5]; ...
                 [-.25 -.5; .25 -.5];[-.25 .5; .25 .5];[-1 .5; 1 .5];[-1 2.5; 1 2.5];     ...
                 [-.25 1.5; .25 1.5];[-1 1.5; 1 1.5]
                };
    coords{2} = {[-.25 -3.5; .25 -3.5];[-.25 -2.5; .25 -2.5];[-2 -.5; 2 -.5];             ...
                 [-3 -.5; 3 -.5];[-2 .5; 2 .5]                                            ...
                 };
    coords{3} = {[-2 3.5; 2 3.5];[-1 3.5; 1 3.5];[-.25 3.5; .25 3.5];[-2 2.5; 2 2.5];     ...
                 [-.25 2.5; .25 2.5];[-2 1.5; 2 1.5];[-2 .5; 2 .5];[-1 -.5; 1 -.5];       ...
                 [-1 -1.5; 1 1.5];[-.25 -1.5; .25 -1.5]                                   ...
                 };
               
    if iscell(filters.mazetype); mazetype = filters.mazetype{1}; else mazetype = filters.mazetype; end
    fn               = sprintf('fullGrid_clusters_%sMaze_%s_%s_perfTh%d',mazetype,filters.epoch,filters.genotype,100*filters.perfTh);
    if isfield(filters,'power'); fn = sprintf('%s_%dmW',fn,round(filters.power)); end
    fn               = sprintf('%s_nBinLogReg%d',fn,cfg.nBins_logReg);

    dataLbl          = sprintf('fullGrid, clusters, maze: %s, epoch: %s, %s mice, perfTh: %d%%',mazetype,filters.epoch,filters.genotype,100*filters.perfTh);
    if isfield(filters,'power'); dataLbl = sprintf('%s, %dmW',dataLbl,round(filters.power)); end

    data             = load('fullGridBilateral.mat','grid');
    grid             = data.grid;

    % loop over locations
    parfor iLoc = 1:numel(coords)
      locIdx        = getMasterLocationIdx(coords{iLoc});
      lsrPerf{iLoc} = analyzeConcatLsrLog(lg,locIdx,filters,cfg,grid);
    end
    lsrPerf         = FDRbatch(lsrPerf);
    lsrPerf         = mouseLsrStats(lsrPerf);

    % plot and save
    rootdir         = ap.getRootDir(spockFlag);
    cd(rootdir)
    save(fn,'lsrPerf','filters','fn','dataLbl','loclbls','coords')

    fprintf('\nANALYZED %s in %1.1f min\n',upper(fn),toc(thist)/60)
    clear grid lsrPerf filters locIdx 

  %% full grid whole trial, 6 mW, vgat, visually guided
  case 'visGuide'
    thist = tic;
    fprintf('-------------------\n\n')

    filters.mazetype                = 'visualGuide';
    filters.perfTh                  = 0.8;
    filters.deleteNonMainMazeTrials = false;
    filters.power                   = 6;
    filters.powerLogic              = '>=';
    filters.deleteNonMainMazeTrials = true;
    filters.epoch                   = 'whole';
    filters.genotype                = 'vgat';
    filters.gridLabel               = {'fullGridBilateral.mat';      ...
                                       'fullLessGridBilateral.mat';  ...     
                                       'posteriorGridBilateral.mat'; ...
                                       'anteriorExtraGridBilateral.mat'};

    % get grid locations, labels
    lg               = concatLogsLsr([],filters,spockFlag);
    if iscell(filters.mazetype); mazetype = filters.mazetype{1}; else mazetype = filters.mazetype; end
    fn               = sprintf('fullGrid_%sMaze_%s_%s_perfTh%d',mazetype,filters.epoch,filters.genotype,100*filters.perfTh);
    if isfield(filters,'power'); fn = sprintf('%s_%dmW',fn,round(filters.power)); end
    fn               = sprintf('%s_nBinLogReg%d',fn,cfg.nBins_logReg);

    dataLbl          = sprintf('fullGrid, maze: %s, epoch: %s, %s mice, perfTh: %d%%',mazetype,filters.epoch,filters.genotype,100*filters.perfTh);
    if isfield(filters,'power'); dataLbl = sprintf('%s, %dmW',dataLbl,round(filters.power)); end

    data             = load('fullGridBilateral.mat','grid');
    grid             = data.grid;
    locIdx           = getMasterLocationIdx(grid);

    % loop over locations
    parfor iLoc = 1:numel(locIdx)
      lsrPerf{iLoc} = analyzeConcatLsrLog(lg,locIdx(iLoc),filters,cfg,grid);
    end
    lsrPerf         = FDRbatch(lsrPerf);
    lsrPerf         = mouseLsrStats(lsrPerf);

    % plot and save
    rootdir         = ap.getRootDir(spockFlag);
    cd(rootdir)
    save(fn,'lsrPerf','filters','fn','dataLbl')

    fprintf('\nANALYZED %s in %1.1f min\n',upper(fn),toc(thist)/60)
    clear filters grid lsrPerf locIdx

  %% full grid whole trial, 6 mW, vgat, memory guided
  case 'memGuide'
    thist = tic;
    fprintf('-------------------\n\n')

    filters.mazetype                = 'memGuide_noTowers';
    filters.perfTh                  = 0.6;
    filters.deleteNonMainMazeTrials = true;
    filters.power                   = 6;
    filters.powerLogic              = '==';
    filters.epoch                   = 'whole';
    filters.genotype                = 'vgat';
    filters.gridLabel               = {'fullGridBilateral.mat';      ...
                                       'fullGridBilateral_higherProb.mat'};

    % get grid locations, labels
    lg               = concatLogsLsr([],filters,spockFlag);
    
    fn               = sprintf('fullGrid_%sMaze_%s_%s_perfTh%d',filters.mazetype,filters.epoch,filters.genotype,100*filters.perfTh);
    if isfield(filters,'power'); fn = sprintf('%s_%dmW',fn,round(filters.power)); end
    fn               = sprintf('%s_nBinLogReg%d',fn,cfg.nBins_logReg);

    dataLbl          = sprintf('fullGrid, maze: %s, epoch: %s, %s mice, perfTh: %d%%',filters.mazetype,filters.epoch,filters.genotype,100*filters.perfTh);
    if isfield(filters,'power'); dataLbl = sprintf('%s, %dmW',dataLbl,round(filters.power)); end

    data             = load('fullGridBilateral.mat','grid');
    grid             = data.grid;
    locIdx           = getMasterLocationIdx(grid);

    % loop over locations
    parfor iLoc = 1:numel(locIdx)
      lsrPerf{iLoc} = analyzeConcatLsrLog(lg,locIdx(iLoc),filters,cfg,grid);
    end

    indicators = {'percCorrect','percCorrect_L','percCorrect_R','bias',       ...
      'speed','dur','percFinish','excessTravel','viewDecode','unusualEvents', ...
      'bias_abs','lapse','viewAngSD'};
    lsrPerf         = FDRbatch(lsrPerf,[],indicators);
    lsrPerf         = mouseLsrStats(lsrPerf);

    % plot and save
    rootdir         = ap.getRootDir(spockFlag);
    cd(rootdir)
    save(fn,'lsrPerf','filters','fn','dataLbl')

    fprintf('\nANALYZED %s in %1.1f min\n',upper(fn),toc(thist)/60)
    clear filters grid lsrPerf locIdx
    
  %% full grid whole trial, 6 mW, wt controls
  case 'WTctrls'
    thist = tic;
    fprintf('-------------------\n\n')

    filters.mazetype                = {'last';'last-1'};
    filters.perfTh                  = 0.6;
    filters.power                   = 6;
    filters.powerLogic              = '>=';
    filters.deleteNonMainMazeTrials = true;
    filters.epoch                   = 'whole';
    filters.genotype                = 'wt';
    filters.gridLabel               = {'fullGridBilateral.mat';      ...
                                       'fullLessGridBilateral.mat';  ...     
                                       'posteriorGridBilateral.mat'; ...
                                       'anteriorExtraGridBilateral.mat'};

    % get grid locations, labels
    lg               = concatLogsLsr([],filters,spockFlag);
    if iscell(filters.mazetype); mazetype = filters.mazetype{1}; else mazetype = filters.mazetype; end
    fn               = sprintf('fullGrid_%sMaze_%s_%s_perfTh%d',mazetype,filters.epoch,filters.genotype,100*filters.perfTh);
    if isfield(filters,'power'); fn = sprintf('%s_%dmW',fn,round(filters.power)); end
    fn               = sprintf('%s_nBinLogReg%d',fn,cfg.nBins_logReg);

    dataLbl          = sprintf('fullGrid, maze: %s, epoch: %s, %s mice, perfTh: %d%%',mazetype,filters.epoch,filters.genotype,100*filters.perfTh);
    if isfield(filters,'power'); dataLbl = sprintf('%s, %dmW',dataLbl,round(filters.power)); end

    data             = load('fullGridBilateral.mat','grid');
    grid             = data.grid;
    locIdx           = getMasterLocationIdx(grid);

    % loop over locations
    parfor iLoc = 1:numel(locIdx)
      lsrPerf{iLoc} = analyzeConcatLsrLog(lg,locIdx(iLoc),filters,cfg,grid);
    end
    lsrPerf         = FDRbatch(lsrPerf);
    lsrPerf         = mouseLsrStats(lsrPerf);

    % plot and save
    rootdir         = ap.getRootDir(spockFlag);
    cd(rootdir)
    save(fn,'lsrPerf','filters','fn','dataLbl')


    fprintf('\nANALYZED %s in %1.1f min\n',upper(fn),toc(thist)/60)
    clear filters grid lsrPerf locIdx

  %% combine reduced grids for epoch-specific
  case 'epochs'
    %% cueHalf1
    thist = tic;
    fprintf('-------------------\n\n')

    filters.mazetype                = {'last';'last-1'};
    filters.perfTh                  = 0.6;
    filters.power                   = 6;
    filters.powerLogic              = '>=';
    filters.deleteNonMainMazeTrials = true;
    filters.epoch                   = 'cueHalf1';
    filters.genotype                = 'vgat';
    filters.gridLabel               = {'veryReducedGridBilateral.mat';    ... 
                                       'veryReducedGridBilateral_v2.mat'; ...
                                       'veryReducedGridBilateral_v3.mat'; ...
                                       'ML10_AP25_bilateralGrid.mat';     ...
                                       'ML10_AP30_bilateralGrid.mat';     ...
                                       'PPCHarveyBilateralGrid.mat';      ...
                                       'ML5_AP-25_bilateralGrid.mat';     ...
                                       };

    % allow for slightly different locations in case of V1, aM2
    loclbls   = {'V1','PPC','RSC','mM2','aM2'};
    coords{1} = {[-2.5 -3.5; 2.5 -3.5];[-3 -3.5; 3 -3.5]};
    coords{2} = {[-1.75 -2; 1.75 -2]};
    coords{3} = {[-.5 -2.5; .5 -2.5]};
    coords{4} = {[-.5 0; .5 0]};
    coords{5} = {[-1 3; 1 3];[-1 2.5; 1 2.5]};
    data      = load('veryReducedGridBilateral_v2.mat','grid');
    grid      = data.grid([4 2 5 3 1]);

    % get grid locations, labels
    if iscell(filters.mazetype); mazetype = filters.mazetype{1}; else mazetype = filters.mazetype; end
    fn               = sprintf('reducedGrid_%sMaze_%s_%s_perfTh%d',mazetype,filters.epoch,filters.genotype,100*filters.perfTh);
    if isfield(filters,'power'); fn = sprintf('%s_%dmW',fn,round(filters.power)); end
    fn               = sprintf('%s_nBinLogReg%d',fn,cfg.nBins_logReg);

    dataLbl          = sprintf('reducedGrid, maze: %s, epoch: %s, %s mice, perfTh: %d%%',mazetype,filters.epoch,filters.genotype,100*filters.perfTh);
    if isfield(filters,'power'); dataLbl = sprintf('%s, %dmW',dataLbl,round(filters.power)); end

    % loop over locations
    for iLoc = 1:numel(coords)
      filt{iLoc}              = filters;
      locIdx{iLoc}            = getMasterLocationIdx(coords{iLoc});
      filt{iLoc}.masterLocIdx = locIdx{iLoc};
    end
    parfor iLoc = 1:numel(coords)
      lg                   = concatLogsLsr([],filt{iLoc},spockFlag);
      lsrPerf{iLoc}        = analyzeConcatLsrLog(lg,locIdx{iLoc},filt{iLoc},cfg,grid);
    end
    lsrPerf         = FDRbatch(lsrPerf);
    lsrPerf         = mouseLsrStats(lsrPerf);

    % plot and save
    rootdir         = ap.getRootDir(spockFlag);
    cd(rootdir)
    save(fn,'lsrPerf','filters','fn','dataLbl','loclbls')

    fprintf('\nANALYZED %s in %1.1f min\n',upper(fn),toc(thist)/60)
    clear lsrPerf filt

    %% cueHalf2
    thist = tic;
    fprintf('-------------------\n\n')

    filters.epoch    = 'cueHalf2';

    % get grid locations, labels
    if iscell(filters.mazetype); mazetype = filters.mazetype{1}; else mazetype = filters.mazetype; end
    fn               = sprintf('reducedGrid_%sMaze_%s_%s_perfTh%d',mazetype,filters.epoch,filters.genotype,100*filters.perfTh);
    if isfield(filters,'power'); fn = sprintf('%s_%dmW',fn,round(filters.power)); end
    fn               = sprintf('%s_nBinLogReg%d',fn,cfg.nBins_logReg);

    dataLbl          = sprintf('reducedGrid, maze: %s, epoch: %s, %s mice, perfTh: %d%%',mazetype,filters.epoch,filters.genotype,100*filters.perfTh);
    if isfield(filters,'power'); dataLbl = sprintf('%s, %dmW',dataLbl,round(filters.power)); end

    % loop over locations
    for iLoc = 1:numel(coords)
      filt{iLoc}              = filters;
      locIdx{iLoc}            = getMasterLocationIdx(coords{iLoc});
      filt{iLoc}.masterLocIdx = locIdx{iLoc};
    end
    parfor iLoc = 1:numel(coords)
      lg                   = concatLogsLsr([],filt{iLoc},spockFlag);
      lsrPerf{iLoc}        = analyzeConcatLsrLog(lg,locIdx{iLoc},filt{iLoc},cfg,grid);
    end
    lsrPerf         = FDRbatch(lsrPerf);
    lsrPerf         = mouseLsrStats(lsrPerf);

    % plot and save
    rootdir         = ap.getRootDir(spockFlag);
    cd(rootdir)
    save(fn,'lsrPerf','filters','fn','dataLbl','loclbls')

    fprintf('\nANALYZED %s in %1.1f min\n',upper(fn),toc(thist)/60)
    clear lsrPerf

    %% delay
    thist = tic;
    fprintf('-------------------\n\n')

    filters.epoch    = 'mem';

    % get grid locations, labels
    if iscell(filters.mazetype); mazetype = filters.mazetype{1}; else mazetype = filters.mazetype; end
    fn               = sprintf('reducedGrid_%sMaze_%s_%s_perfTh%d',mazetype,filters.epoch,filters.genotype,100*filters.perfTh);
    if isfield(filters,'power'); fn = sprintf('%s_%dmW',fn,round(filters.power)); end
    fn               = sprintf('%s_nBinLogReg%d',fn,cfg.nBins_logReg);

    dataLbl          = sprintf('reducedGrid, maze: %s, epoch: %s, %s mice, perfTh: %d%%',mazetype,filters.epoch,filters.genotype,100*filters.perfTh);
    if isfield(filters,'power'); dataLbl = sprintf('%s, %mW',dataLbl,round(filters.power)); end

    % loop over locations
    for iLoc = 1:numel(coords)
      filt{iLoc}              = filters;
      locIdx{iLoc}            = getMasterLocationIdx(coords{iLoc});
      filt{iLoc}.masterLocIdx = locIdx{iLoc};
    end
    parfor iLoc = 1:numel(coords)
      lg                   = concatLogsLsr([],filt{iLoc},spockFlag);
      lsrPerf{iLoc}        = analyzeConcatLsrLog(lg,locIdx{iLoc},filt{iLoc},cfg,grid);
    end
    lsrPerf         = FDRbatch(lsrPerf);
    lsrPerf         = mouseLsrStats(lsrPerf);

    % plot and save
    rootdir         = ap.getRootDir(spockFlag);
    cd(rootdir)
    save(fn,'lsrPerf','filters','fn','dataLbl','loclbls')
    
    fprintf('\nANALYZED %s in %1.1f min\n',upper(fn),toc(thist)/60)
    clear lsrPerf filters filt locIdx coords grid

end

%% shut down parallel pool
delete(poolobj);
catch ME
  displayException(ME)
end
