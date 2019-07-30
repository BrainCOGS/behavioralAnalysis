function motorErr = analyzeMotorErrors(lg,cfg)

% motorErr = analyzeMotorErrors(lg,cfg)
% incidence of motor errors from flattened log, lg. cfg is data structure
% with analysis params, optional

%%
if nargin < 2
  cfg.minNumTrials = 1000;
  cfg.motorEvents  = {'cueTurnAround','earlyTurn','enteredOppositeArm','lowSpeed','excessTravel','any'};
  cfg.deltaTh      = 5;
end

if sum(strcmpi(cfg.motorEvents,'any')) == 0
  cfg.motorEvents{end+1} = 'any';
end

% just get mice with minNT, other criteria shouldn' apply
[~, info]        = cleanupConcatLog(lg, cfg.minNumTrials);
mouseID          = info.mouseIndices;
nmice            = numel(mouseID);
nevents          = numel(cfg.motorEvents);

%%
motorErr.eventLbls         = cfg.motorEvents;
motorErr.deltaTh           = cfg.deltaTh;
motorErr.mouseID           = info.mouseIndices;
motorErr.freqOverall       = nan(nmice,nevents);
motorErr.freqEasyDelta     = nan(nmice,nevents);
motorErr.freqHardDelta     = nan(nmice,nevents);
motorErr.freqCorrectTrials = nan(nmice,nevents);
motorErr.freqErrorTrials   = nan(nmice,nevents);

%%
for iMouse = 1:nmice
  
  trialID = cell(nevents,1);
  for iEvent = 1:nevents
    [trialID{iEvent}, delta, choice]          = unusualMotorEvents(lg, mouseID(iMouse), cfg.motorEvents{iEvent});
    motorErr.freqOverall(iMouse,iEvent)       = 100 * sum(trialID{iEvent})/numel(trialID{iEvent});
    motorErr.freqEasyDelta(iMouse,iEvent)     = 100 * sum(trialID{iEvent} & abs(delta) > cfg.deltaTh)/sum(abs(delta) > cfg.deltaTh);
    motorErr.freqHardDelta(iMouse,iEvent)     = 100 * sum(trialID{iEvent} & abs(delta) < cfg.deltaTh)/sum(abs(delta) > cfg.deltaTh);
    correctTrials                             = (delta > 0 & choice == analysisParams.rightCode) | (delta < 0 & choice == analysisParams.leftCode);
    errorTrials                               = (delta > 0 & choice == analysisParams.leftCode)  | (delta < 0 & choice == analysisParams.rightCode);
    motorErr.freqCorrectTrials(iMouse,iEvent) = 100 * sum(trialID{iEvent} & correctTrials)/sum(correctTrials);
    motorErr.freqErrorTrials(iMouse,iEvent)   = 100 * sum(trialID{iEvent} & errorTrials)/sum(errorTrials);
  end
  
  % overall rate
  anyevent                              = strcmpi(cfg.motorEvents,'any');
  speedid                               = strcmpi(cfg.motorEvents,'lowSpeed');
  motorErr.overallRate_allTypes(iMouse) = 100 .* sum(trialID{anyevent==1} & ~trialID{speedid==1})/numel(trialID{anyevent==1}); 
  
end

%% stats
% 2-way ANOVA with repeated measures
events_hard        = motorErr.freqHardDelta(:,~anyevent);
events_easy        = motorErr.freqEasyDelta(:,~anyevent);
events_correct     = motorErr.freqCorrectTrials(:,~anyevent);
events_error       = motorErr.freqErrorTrials(:,~anyevent);
[nmice, nevents]   = size(events_hard);

data_diff          = [];
data_outc          = [];
mousevec           = [];
eventvec           = [];
diffvec            = [];
outvec             = [];
for iEvent = 1:nevents
  data_diff(end+1:end+nmice,:)  = events_easy(:,iEvent);
  data_diff(end+1:end+nmice,:)  = events_hard(:,iEvent);
  data_outc(end+1:end+nmice,:)  = events_correct(:,iEvent);
  data_outc(end+1:end+nmice,:)  = events_error(:,iEvent);
  mousevec(end+1:end+nmice*2,:) = repmat((1:nmice)',[2 1]);
  eventvec(end+1:end+nmice*2,:) = repmat(iEvent.*ones(nmice,1),[2 1]);
  diffvec(end+1:end+nmice*2,:)  = [ones(nmice,1); 2.*ones(nmice,1)];
  outvec(end+1:end+nmice*2,:)   = [ones(nmice,1); 2.*ones(nmice,1)];
end

motorErr.stats.testname               = '2-way RM ANOVA';
motorErr.stats.varname_difficulty     = {'eventType','difficulty','mouseid'};
motorErr.stats.varname_errVScorrect   = {'eventType','outcome','mouseid'};

[motorErr.stats.pvals_anova_diffic,motorErr.stats.table_anova_diffic,st]        ...
                                      =                                         ...
                           anovan(data_diff,{eventvec, diffvec, mousevec},      ...
                                  'varnames',motorErr.stats.varname_difficulty,'display','off');
motorErr.stats.multCompTable_diffic   = multcompare(st,'display','off','dimension',[1 2]);
motorErr.stats.multCompPvals_diffic   = motorErr.stats.multCompTable_diffic(:,end);

[motorErr.stats.pvals_anova_err,motorErr.stats.table_anova_err,st]              ...
                                      =                                         ...
                           anovan(data_outc,{eventvec, outvec, mousevec},       ...
                                  'varnames',motorErr.stats.varname_errVScorrect,'display','off');
motorErr.stats.multCompTable_err      = multcompare(st,'display','off');
motorErr.stats.multCompPvals_err      = motorErr.stats.multCompTable_err(:,end);

%% sign rank tests
for iEvent = 1:nevents
  motorErr.stats.easyVsHard_pval(iEvent)    = signrank(motorErr.freqEasyDelta(:,iEvent),     ...
                                                       motorErr.freqHardDelta(:,iEvent));
  motorErr.stats.correctVsErr_pval(iEvent)  = signrank(motorErr.freqCorrectTrials(:,iEvent), ...
                                                       motorErr.freqErrorTrials(:,iEvent));
end

%% pairwise comparisons within events using multi comp
for iEvent = 1:nevents
  rowidx = motorErr.stats.multCompTable_diffic(:,1) == iEvent & ...
           motorErr.stats.multCompTable_diffic(:,2) == iEvent+nevents;
  motorErr.stats.easyVsHard_pval_multComp(iEvent)    = motorErr.stats.multCompPvals_err(rowidx);
  motorErr.stats.correctVsErr_pval_multComp(iEvent)  = motorErr.stats.multCompPvals_diffic(rowidx);
end
