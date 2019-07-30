function taskComp = compareTasksMotor(lg,lg_visGuide,lg_memGuide,analysisFilePath,plotFlag,justWFmice)

% taskComp = compareTasksMotor(lg,lg_visGuide,lg_memGuide,analysisFilePath,plotFlag,justWFmice)
% compares motor indicators across different tasks

%% initalize
if nargin < 3; lg_memGuide      = [];                      end
if nargin < 4; analysisFilePath = analysisParams.savepath; end
if nargin < 5; plotFlag         = false;                   end
if nargin < 6; justWFmice       = false;                   end

cfg.visGuideExcludeID   = nan;%22; % this mouse was not used in towers
cfg.towersCl            = widefieldParams.darkgray; %[0 0 0];
cfg.ctrlCl              = widefieldParams.darkgreen; %[.6 .6 .6];
cfg.posBins             = 0:5:300;
cfg.posBinsVa           = 0:1:300;
cfg.posBinRange         = [-20, 50];
cfg.numSides            = 2;
cfg.cueVisAt            = 10;

%% get combined log for baseline behavior (opto & wf)
lg.firstTrialofBlock = [true diff(lg.meanPerfBlock)~=0];
lg_dirty             = lg;
lg                   = cleanupConcatLog(lg, 0);

%%
lg_visGuide.firstTrialofBlock = [true diff(lg_visGuide.meanPerfBlock)~=0];
lg_visGuide_dirty             = lg_visGuide;
lg_visGuide                   = cleanupConcatLog(lg_visGuide, 0, 0, 0, 0, 0, .8);

%%
if ~isempty(lg_memGuide)
  lg_memGuide.firstTrialofBlock = [true diff(lg_memGuide.meanPerfBlock)~=0];
  lg_memGuide_dirty             = lg_memGuide;
  lg_memGuide                   = cleanupConcatLog(lg_memGuide, 0);
end
%% cleanup
if justWFmice
  mid = unique(lg_visGuide.mouseID);
  mid = mid(end-5:end);
  idx = ~ismember(lg_visGuide.mouseID,mid);
  fieldls = fields(lg_visGuide);
  for iF = 1:numel(fieldls)
    if strcmpi(fieldls{iF},'mouseNames'); continue; end
    lg_visGuide.(fieldls{iF})(idx)      = [];
  end
  
  idx = ~ismember(lg.mouseID,mid);
  fieldls = fields(lg);
  for iF = 1:numel(fieldls)
    if strcmpi(fieldls{iF},'mouseNames'); continue; end
    lg.(fieldls{iF})(idx)      = [];
  end
else
  fieldls = fields(lg_visGuide);
  idx     = lg_visGuide.mouseID == cfg.visGuideExcludeID;
  for iF = 1:numel(fieldls)
    if strcmpi(fieldls{iF},'mouseNames'); continue; end
    lg_visGuide.(fieldls{iF})(idx)      = [];
  end
end

%% view angle SD
thisVSD                       = xMouseViewAngSD(lg,0,cfg.posBins);
taskComp.towers.viewAngSD     = thisVSD.viewAngSD;
thisVSD                       = xMouseViewAngSD(lg_visGuide,0,cfg.posBins);
taskComp.visGuide.viewAngSD   = thisVSD.viewAngSD;
if ~isempty(lg_memGuide)
  thisVSD                     = xMouseViewAngSD(lg_memGuide,0,cfg.posBins);
  taskComp.memGuide.viewAngSD = thisVSD.viewAngSD;
end

%% motor errors
taskComp.towers.motorErr      = xMouseMotorErrors(lg_dirty);
taskComp.visGuide.motorErr    = xMouseMotorErrors(lg_visGuide_dirty);
if ~isempty(lg_memGuide)
  taskComp.memGuide.motorErr  = xMouseMotorErrors(lg_memGuide_dirty);
end
%% speed (overall and by section)
taskComp.towers.speed         = xMouseSpeed(lg,'all');
taskComp.visGuide.speed       = xMouseSpeed(lg_visGuide,'all');
taskComp.towers.speed_cue     = xMouseSpeed(lg,'cue');
taskComp.visGuide.speed_cue   = xMouseSpeed(lg_visGuide,'cue');
taskComp.towers.speed_del     = xMouseSpeed(lg,'delay');
taskComp.visGuide.speed_del   = xMouseSpeed(lg_visGuide,'delay');

if ~isempty(lg_memGuide)
  taskComp.memGuide.speed       = xMouseSpeed(lg_memGuide,'all');
  taskComp.memGuide.speed_cue   = xMouseSpeed(lg_memGuide,'cue');
  taskComp.memGuide.speed_del   = xMouseSpeed(lg_memGuide,'delay');
end

%% acceleration (by section)
taskComp.towers.accel_first50   = xMouseAcceleration(lg,[0 50]);
taskComp.visGuide.accel_first50 = xMouseAcceleration(lg_visGuide,[0 50]);
taskComp.towers.accel_last50    = xMouseAcceleration(lg,[250 300]);
taskComp.visGuide.accel_last50  = xMouseAcceleration(lg_visGuide,[250 300]);

if ~isempty(lg_memGuide)
  taskComp.memGuide.accel_first50 = xMouseAcceleration(lg_memGuide,[0 50]);
  taskComp.memGuide.accel_last50  = xMouseAcceleration(lg_memGuide,[250 300]);
end

%% turn onset
taskComp.towers.turnOnset       = xMouseTurnOnset(lg);
taskComp.visGuide.turnOnset     = xMouseTurnOnset(lg_visGuide);

if ~isempty(lg_memGuide)
  taskComp.memGuide.turnOnset   = xMouseTurnOnset(lg_memGuide);
end

%% plot and do paired stats
vecs    = fieldnames(taskComp.towers);

if plotFlag
  ap      = analysisParams;
  [nr,nc] = subplotOrg(numel(vecs),4);
  fh      = figure;
  ap.applyFigDefaults(fh,[nc nr],'w');
end

if isempty(lg_memGuide)
  for iVec = 1:numel(vecs)
    thisAT = taskComp.towers.(vecs{iVec});
    thisVG = taskComp.visGuide.(vecs{iVec});
    taskComp.towers.([vecs{iVec} '_avg'])   = nanmean(thisAT);
    taskComp.towers.([vecs{iVec} '_sem'])   = nanstd(thisAT)./sqrt(numel(thisAT)-1);
    taskComp.visGuide.([vecs{iVec} '_avg']) = nanmean(thisVG);
    taskComp.visGuide.([vecs{iVec} '_sem']) = nanstd(thisVG)./sqrt(numel(thisVG)-1);

    if size(thisAT,1) < size(thisAT,2)
      thisAT = thisAT'; 
      thisVG = thisVG';
    end
    if lillietest(thisAT) || lillietest(thisVG)
      if jusWFmice
        taskComp.stats.(['p_' vecs{iVec}])     = signrank(thisAT,thisVG);
        taskComp.stats.(['test_' vecs{iVec}])  = 'signrank';
      else
        taskComp.stats.(['p_' vecs{iVec}])     = ranksum(thisAT,thisVG);
        taskComp.stats.(['test_' vecs{iVec}])  = 'ranksum';
      end
    else
      if justWFmice
        [~,taskComp.stats.(['p_' vecs{iVec}])] = ttest(thisAT,thisVG);
        taskComp.stats.(['test_' vecs{iVec}])  = 'ttest';
      else
        [~,taskComp.stats.(['p_' vecs{iVec}])] = ttest2(thisAT,thisVG);
        taskComp.stats.(['test_' vecs{iVec}])  = 'ttest2';
      end
    end

    if plotFlag
      subplot(nr,nc,iVec); hold on
      thisAT = thisAT';
      thisVG = thisVG';
      plot([ones(size(thisAT)); 1+ones(size(thisAT))],[thisAT; thisVG], ...
         '-','linewidth',.5,'color',[.6 .6 .6])
      errorbar(0.85,nanmean(thisAT),nanstd(thisAT)./sqrt(numel(thisAT)-1), ...
           '.-','linewidth',.5,'color',cfg.towersCl,'markersize',20) 
      errorbar(2.15,nanmean(thisVG),nanstd(thisVG)./sqrt(numel(thisVG)-1), ...
           '.-','linewidth',.5,'color',cfg.ctrlCl,'markersize',20) 

      yl = get(gca,'ylim');
      text(1.5,yl(2)*.95,sprintf('P = %1.2g',taskComp.stats.(['p_' vecs{iVec}])), ...
           'color','k','horizontalAlignment','center','fontsize',8)

      set(gca,'xtick',1:2,'xticklabel',{'accum. towers';'control'})
      rotateXLabels(gca,45)
      ylabel(removeUnderscores(vecs{iVec}),'fontsize',12)
      xlim([.65 2.35])
      box off
    end
  end
else
  for iVec = 1:numel(vecs)
    thisAT = taskComp.towers.(vecs{iVec});
    thisVG = taskComp.visGuide.(vecs{iVec});
    thisMG = taskComp.memGuide.(vecs{iVec});
    taskComp.towers.([vecs{iVec} '_avg'])   = nanmean(thisAT);
    taskComp.towers.([vecs{iVec} '_sem'])   = nanstd(thisAT)./sqrt(numel(thisAT)-1);
    taskComp.visGuide.([vecs{iVec} '_avg']) = nanmean(thisVG);
    taskComp.visGuide.([vecs{iVec} '_sem']) = nanstd(thisVG)./sqrt(numel(thisVG)-1);
    taskComp.memGuide.([vecs{iVec} '_avg']) = nanmean(thisMG);
    taskComp.memGuide.([vecs{iVec} '_sem']) = nanstd(thisMG)./sqrt(numel(thisVG)-1);

    if size(thisAT,1) < size(thisAT,2)
      thisAT = thisAT'; 
      thisVG = thisVG';
      thisMG = thisMG';
    end
    
    taskComp.stats.(['p_' vecs{iVec}])    = ...
      anovan([thisVG;thisMG;thisAT],{[ones(size(thisVG)); ones(size(thisMG))+1; ones(size(thisAT))+2]},'display','off');
    taskComp.stats.(['test_' vecs{iVec}]) = 'anova';
  end
end

%% cue-triggered view angle
taskComp.towers.towerTrigVAng   = cueTrigViewAngle(lg,cfg);  
taskComp.visGuide.towerTrigVAng = cueTrigViewAngle(lg_visGuide,cfg); 

% 2-way RM ANOVA for task and spatial view angle modulation
idx          = taskComp.towers.towerTrigVAng.posAxis > 0;
data_vg      = taskComp.visGuide.towerTrigVAng.alignedVAng_sideAvg_mouse(idx,:);
data_tow     = taskComp.towers.towerTrigVAng.alignedVAng_sideAvg_mouse(idx,:);
[npos,nmice] = size(data_vg);
[npost,nmicet] = size(data_tow);
% mousevec     = repmat((1:nmice),[npos 1]);
% mousevec     = repmat(mousevec(:),[2 1]);
posvec       = repmat((1:npos)',[1 nmice]);
posvect      = repmat((1:npost)',[1 nmicet]);
posvec       = [posvec(:); posvect(:)];
taskvec      = [ones(npos*nmice,1); 1+ones(npost*nmicet,1)];
datavec      = [data_vg(:); data_tow(:)];

[taskComp.stats.towerTrigVAng_ANOVA_p,taskComp.stats.towerTrigVAng_ANOVA_table,taskComp.stats.towerTrigVAng_ANOVA_stats] ...
             = anovan(datavec,{taskvec,posvec},'varnames',{'task','position'},'display','off');
           
%% save
if ~isempty(analysisFilePath)
  cd(analysisFilePath)
  save taskCompMotor taskComp
  
  if plotFlag 
    saveas(gcf,'taskCompMotor')
  end

end

end

%% -----------------------------------------
%% tower-triggered view angles
function va = cueTrigViewAngle(lg,cfg)

% Sample view angle at various positions
viewAngle               = sampleViewAngleVsY(lg.pos, cfg.posBinsVa);
mid                     = unique(lg.mouseID);

% Align and average view angle to left/right cue onsets
nPosBins                = diff(cfg.posBinRange) + 1;
relPosBins              = cfg.posBinRange(1):cfg.posBinRange(end);
alignedVAng             = nan(nPosBins, cfg.numSides, numel(mid));    % (pos,side,mouse)
for iMouse = 1:numel(mid)
  %% Trial selection 
  [choice, cuePos_L, cuePos_R, sel]       ...
                        = selectMouseTrials(lg, mid(iMouse), 'choice', 'cuePos_L', 'cuePos_R');
  vAngle                = viewAngle(:,sel);
  cueOnset              = {cuePos_L{:}; cuePos_R{:}};

  %% Subtract mean trajectory for a given choice
  for iChoice = 1:cfg.numSides
    sel                 = choice == (iChoice - 1);
    meanVAng            = mean(vAngle(:,sel), 2);
    vAngle(:,sel)       = bsxfun(@minus, vAngle(:,sel), meanVAng);
  end

  % Average view angle triggered to cues on a particular side
  for iSide = 1:cfg.numSides
    cueTrigVAng         = nan(nPosBins, numel(choice));
    for iTrial = 1:numel(choice)
      cues              = cueOnset{iSide,iTrial} - cfg.cueVisAt;
      iCuePos           = binarySearch(cfg.posBinsVa, cues, 0, 2);
      for iCue = 1:numel(cues)
        % Range of position bins aligned to the cue location
        srcRange        = iCuePos(iCue) + relPosBins;
        if srcRange(1) < 1
          iFirst        = -srcRange(1) + 2;
        else
          iFirst        = 1;
        end
        if srcRange(end) > size(vAngle,1)
          iLast         = numel(srcRange) - (size(vAngle,1) - srcRange(end)) - 1;
        else
          iLast         = numel(srcRange);
        end

        % Copy view angle within valid range
        tgtRange        = iFirst:iLast;
        cueTrigVAng(tgtRange,iTrial)  = vAngle(srcRange(tgtRange),iTrial);
      end
    end

    %% Record for this side of cue and this mouse
    alignedVAng(:,iSide,iMouse)       = mean(cueTrigVAng, 2, 'omitnan');
  end
end

sideAvg                      = alignedVAng;
sideAvg(:,2,:)               = -sideAvg(:,2,:); % make it such that ipsiversive is positive
sideAvg                      = nanmean(sideAvg,2); % avg left and right
basel                        = nanmean(sideAvg(relPosBins<0,:,:));
sideAvg                      = sideAvg - repmat(basel,[size(sideAvg,1) 1 1]);
va.alignedVAng_mouse         = alignedVAng;
va.alignedVAng_avg           = nanmean(alignedVAng,3);
va.alignedVAng_sem           = nanstd(alignedVAng,0,3)./sqrt(numel(mid)-1);
va.alignedVAng_sideAvg_mouse = sideAvg;
va.alignedVAng_sideAvg_avg   = squeeze(nanmean(sideAvg,3));
va.alignedVAng_sideAvg_sem   = squeeze(nanstd(sideAvg,0,3))./sqrt(numel(mid)-1);
va.posAxis                   = relPosBins;

end