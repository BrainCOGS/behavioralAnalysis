function concatLogs_subtrial

%%
spockFlag                = isThisSpock;
ap                       = analysisParams;
rootdir                  = ap.getRootDir(spockFlag);

%% cfg
cfg.areaOrder           = {'V1','mV2','PPC','RSC','Post','mM2','aM2','M1','Front'};
cfg.epochOrder          = {'cueQuart1','cueQuart2','cueHalf1','cueQuart3','cueHalf2','mem'};
cfg.minNumTrialsLogReg  = 500;
cfg.nBins               = 4;
cfg.minNumTrials        = 10;

%% filters
filters.epoch                   = {'cueHalf1','cueHalf2','mem','cueQuart1','cueQuart2','cueQuart3'};
filters.mazetype                = {'last';'last-1'};
filters.perfTh                  = 0.6;
filters.power                   = 6;
filters.powerLogic              = '<=';
filters.deleteNonMainMazeTrials = true;
filters.genotype                = 'vgat';
filters.gridLabel               = {'veryReducedGridBilateral.mat';          ...
                                   'veryReducedGridBilateral_v2.mat';       ...
                                   'veryReducedGridBilateral_v3.mat';       ...
                                   'ML10_AP25_bilateralGrid.mat';           ...
                                   'ML10_AP30_bilateralGrid.mat';           ...
                                   'PPCHarveyBilateralGrid.mat';            ...
                                   'ML5_AP-25_bilateralGrid.mat';           ...
                                   'veryReducedGridBilateral_mM2V1mV2.mat'; ...
                                   'veryReducedGridBilateral_v2.mat';       ...
                                   'veryReducedGridBilateral_v4.mat';       ...
                                   'veryReducedGridBilateral_v5.mat';       ...
                                   'multiAreaGridBilateral.mat';            ...
                                   'veryReducedGrid.mat';                   ...
                                   'veryReducedGridUnilateral_v2.mat';      ...
                                   'veryReducedGridUnilateral_v5.mat';      ...
                                   'fullGridBilateral.mat';                 ...
                                   'fullLessGridBilateral.mat';             ...
                                   'posteriorGridBilateral.mat';            ...
                                   'anteriorExtraGridBilateral.mat';        ...
                                   'mV2_M1_Post_Front_bilateralGrid.mat';   ...
                                   'mV2_mM2_bilateralGrid.mat'
                                   };

% bilateral and unilateral locations
loclbls   = {'V1','mV2','PPC','RSC','mM2','aM2','M1','Post','Front'};

% bilateral
coords{1}     = [-3 -3.5; 3 -3.5];
coords{end+1} = [-2 -2.5; 2 -2.5];
coords{end+1} = [-1.75 -2; 1.75 -2];
coords{end+1} = [-.5 -2.5; .5 -2.5];
coords{end+1} = [-.5 0; .5 0];
coords{end+1} = [-1 3; 1 3];
coords{end+1} = [-2 1; 2 1];
coords{end+1} = [-3 -3.5; -1.75 -2; -0.5 -2.5; 0.5 -2.5; 1.75 -2; 3 -3.5];
coords{end+1} = [-1 3; -2 1; -0.5 0; 0.5 0; 2 1; 1 3];

locIdx               = getMasterLocationIdx(coords');
% filters.masterLocIdx = locIdx;

%% get full log
lg                   = concatLogsLsr([],filters,spockFlag);

%% add location labels, epoch ranges
lg.loclbl               = cell(1,numel(lg.choice));
lg.laserEpoch_yposRange = cell(1,numel(lg.choice));
for iTrial = 1:numel(lg.choice)
  if ~lg.laserON(iTrial)
    lg.loclbl{iTrial}               = 'none';
    lg.laserEpoch_yposRange{iTrial} = single([nan nan]);
  else
    % loc label
    subidx     = find(cellfun(@(x)(numel(x)),locIdx) == numel(lg.galvoPosIdxMaster{iTrial}));
    thisLocIdx = cellfun(@(x)(sum(x==lg.galvoPosIdxMaster{iTrial})==numel(lg.galvoPosIdxMaster{iTrial})),locIdx(subidx));
    if sum(thisLocIdx) == 0
      lg.loclbl{iTrial} = 'unknown';
    else
      lg.loclbl{iTrial} = loclbls{subidx(thisLocIdx)};
    end

    % epoch range
    switch lg.laserEpoch{iTrial}
      case 'cueHalf1'
        lg.laserEpoch_yposRange{iTrial} = single([0 100]);
      case 'cueHalf2'
        lg.laserEpoch_yposRange{iTrial} = single([100 200]);
      case 'mem'
        lg.laserEpoch_yposRange{iTrial} = single([200 300]);
      case 'cueQuart1'
        lg.laserEpoch_yposRange{iTrial} = single([0 50]);
      case 'cueQuart2'
        lg.laserEpoch_yposRange{iTrial} = single([50 100]);
      case 'cueQuart3'
        lg.laserEpoch_yposRange{iTrial} = single([100 150]);
      case 'memHalf1'
        lg.laserEpoch_yposRange{iTrial} = single([200 250]);
      case 'cue'
        lg.laserEpoch_yposRange{iTrial} = single([0 200]);
      case 'whole'
        lg.laserEpoch_yposRange{iTrial} = single([0 300]);
    end
  end
end

%% delete undesired trials / blocks
idx         = lg.meanPerfBlockCtrl < filters.perfTh | strcmpi(lg.loclbl{iTrial},'unknown') | ...
              strcmpi(lg.laserEpoch{iTrial},'3epochs') | strcmpi(lg.laserEpoch{iTrial},'memHalf1') | ...
              strcmpi(lg.laserEpoch{iTrial},'cue') | strcmpi(lg.laserEpoch{iTrial},'whole');

vectorNames = fieldnames(lg);
for iVec = 1:numel(vectorNames)
  if isfield(lg,vectorNames{iVec}) && ~strcmpi(vectorNames{iVec},'keyFrameLabels')
    lg.(vectorNames{iVec})(idx) = [];
  end
end

lg.firstTrialofBlock = [true diff(lg.meanPerfBlock)~=0];
lg                   = cleanupConcatLog(lg, 0);

% ---
%% now create new log with some metainfo and processed view angles
% copy only these fields
newlg.rightCode            = 1;
newlg.leftCode             = 0;
newlg.ypos_vals            = 0:300;
newlg.epochList            = epochList;
newlg.inactLocationList    = loclbls;
newlg.choice               = lg.choice;
newlg.trialType            = lg.trialType;
newlg.cuePos_R             = lg.cuePos_R;
newlg.cuePos_L             = lg.cuePos_L;
newlg.nCues_RminusL        = cellfun(@numel,lg.cuePos_R) - cellfun(@numel,lg.cuePos_L);
newlg.mouseID              = lg.mouseID;
newlg.sessionID            = lg.sessionID;
newlg.laserON              = lg.laserON;
newlg.laserEpoch           = lg.laserEpoch;
newlg.loclbl               = lg.loclbl;
newlg.laserEpoch_yposRange = lg.laserEpoch_yposRange;


%% downsample and mean-subtract view angle
viewAngle           = sampleViewAngleVsY(lg.pos, newlg.ypos_vals);

%%
viewAngle_meanSubt  = single(zeros(size(viewAngle)));
mice                = unique(lg.mouseID);
for iMouse = 1:numel(mice)
  idx = lg.mouseID == mice(iMouse);
  viewAngle_meanSubt(:,idx) = viewAngle(:,idx) - nanmean(viewAngle(:,idx),2);
end

%%
newlg.viewAngle_byYpos = cell(1,numel(newlg.choice));
for iTrial = 1:numel(newlg.viewAngle_byYpos)
  newlg.viewAngle_byYpos{iTrial} = viewAngle_meanSubt(:,iTrial);
end

%% cleanup large view angle and deleted towers trials
isCtrlTrial = (cellfun(@(x,y)(sum(x > 100) == 0 & sum(y > 100) == 0),lg.cuePos_R,lg.cuePos_L) | ...
              cellfun(@(x,y)(sum(x < 100) == 0 & sum(y < 100) == 0),lg.cuePos_R,lg.cuePos_L)) & ...
              (lg.nCues_RminusL == 0 | (lg.nCues_RminusL > 0 & lg.trialType == analysisParams.leftCode) | ...
              (lg.nCues_RminusL < 0 & lg.trialType == analysisParams.rightCode)) & ~lg.laserON;
largeVa     = cellfun(@(x)(any(abs(x(1:250))>60) | any(abs(x(1:200))>40)),newlg.viewAngle_byYpos);
badtrial    = largeVa | isCtrlTrial;

fieldls     = fields(newlg);
for iField = 1:numel(fieldls)
  if strcmp(fieldls{iField},'rightCode') || strcmp(fieldls{iField},'leftCode') || ...
     strcmp(fieldls{iField},'ypos_vals') || strcmp(fieldls{iField},'epochList') || strcmp(fieldls{iField},'inactLocationList')
   continue
  end
  newlg.(fieldls{iField})(badtrial) = [];
end

%% save
lg = newlg;
cd(rootdir)
save concatLog_subtrial_inactivation lg
