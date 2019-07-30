function th = analyzeTrialHistory(logfp, cfg)

% th = analyzeTrialHistory(logfp, cfg)
% analyzes choice and rewward history influences on behavior
% logfp can be either the flattened log or its path
% cfg is optional analysis parameters structure

if nargin < 1
  logfp = [analysisParams.savepathBehav 'concatLog_lastMaze_perTh60_noLaserSess.mat'];
end
if nargin < 2
  cfg.minNumTrials = 1000;
  cfg.nback        = 5;
  cfg.nConsecPers  = 3;
end

if isstruct(logfp)
  lg = logfp;
else
  load(logfp,'lg')
end

%% clean up log, keeping history trials
[lg,info] = cleanupConcatLog(lg, cfg.minNumTrials, false, cfg.nback, true, false);
th.info   = info;
th.nback  = cfg.nback;

%% ------------------------------------------------------------------------
%% METAMOUSE
%% ------------------------------------------------------------------------

%% vectors
choice      = lg.choice;
delta       = lg.nCues_RminusL;
trialtype   = lg.trialType;
mouseid     = lg.mouseID;
reward      = choice == trialtype;
error       = choice ~= trialtype;
wentR       = choice == analysisParams.rightCode;
wentL       = choice == analysisParams.leftCode;

%% 1 trial back vectors
prev_wR     = [0 wentR(1:end-1)];
prev_wL     = [0 wentL(1:end-1)];
prev_rw     = [0 reward(1:end-1)];
prev_err    = [0 error(1:end-1)];
prev_delta  = [0 delta(1:end-1)];

prev_hardR  = prev_delta > 0 & prev_delta <= 5;
prev_easyR  = prev_delta > 5;
prev_hardL  = prev_delta < 0 & prev_delta >= -5;
prev_easyL  = prev_delta < -5;

%% one trial back, action + outcome
th.psych   = psychometricFit(choice,delta,0,[],0);
th.psychC  = psychometricFit(choice(prev_rw==1),delta(prev_rw==1),0,[],0);
th.psychE  = psychometricFit(choice(prev_err==1),delta(prev_err==1),0,[],0);
th.psychRC = psychometricFit(choice(prev_rw & prev_wR),delta(prev_rw & prev_wR),0,[],0);
th.psychLC = psychometricFit(choice(prev_rw & prev_wL),delta(prev_rw & prev_wL),0,[],0);
th.psychRE = psychometricFit(choice(prev_err & prev_wR),delta(prev_err & prev_wR),0,[],0);
th.psychLE = psychometricFit(choice(prev_err & prev_wL),delta(prev_err & prev_wL),0,[],0);

%% one trial back, just delta
th.psychEasyR       = psychometricFit(choice(prev_easyR),delta(prev_easyR),0,[],0);
th.psychHardR       = psychometricFit(choice(prev_hardR),delta(prev_hardR),0,[],0);
th.psychEasyL       = psychometricFit(choice(prev_easyL),delta(prev_easyL),0,[],0);
th.psychHardL       = psychometricFit(choice(prev_hardL),delta(prev_hardL),0,[],0);

%% one trial back, action + delta
th.psychEasyRwentR  = psychometricFit(choice(prev_easyR & prev_wR),delta(prev_easyR & prev_wR),0,[],0);
th.psychHardRwentR  = psychometricFit(choice(prev_hardR & prev_wR),delta(prev_hardR & prev_wR),0,[],0);
th.psychEasyLwentR  = psychometricFit(choice(prev_easyL & prev_wR),delta(prev_easyL & prev_wR),0,[],0);
th.psychHardLwentR  = psychometricFit(choice(prev_hardL & prev_wR),delta(prev_hardL & prev_wR),0,[],0);

th.psychEasyRwentL  = psychometricFit(choice(prev_easyR & prev_wL),delta(prev_easyR & prev_wL),0,[],0);
th.psychHardRwentL  = psychometricFit(choice(prev_hardR & prev_wL),delta(prev_hardR & prev_wL),0,[],0);
th.psychEasyLwentL  = psychometricFit(choice(prev_easyL & prev_wL),delta(prev_easyL & prev_wL),0,[],0);
th.psychHardLwentL  = psychometricFit(choice(prev_hardL & prev_wL),delta(prev_hardL & prev_wL),0,[],0);

%% one trial back, outcome + delta
th.psychEasyRrew    = psychometricFit(choice(prev_easyR & prev_rw),delta(prev_easyR & prev_rw),0,[],0);
th.psychHardRrew    = psychometricFit(choice(prev_hardR & prev_rw),delta(prev_hardR & prev_rw),0,[],0);
th.psychEasyLrew    = psychometricFit(choice(prev_easyL & prev_rw),delta(prev_easyL & prev_rw),0,[],0);
th.psychHardLrew    = psychometricFit(choice(prev_hardL & prev_rw),delta(prev_hardL & prev_rw),0,[],0);

th.psychEasyRerr    = psychometricFit(choice(prev_easyR & prev_err),delta(prev_easyR & prev_err),0,[],0);
th.psychHardRerr    = psychometricFit(choice(prev_hardR & prev_err),delta(prev_hardR & prev_err),0,[],0);
th.psychEasyLerr    = psychometricFit(choice(prev_easyL & prev_err),delta(prev_easyL & prev_err),0,[],0);
th.psychHardLerr    = psychometricFit(choice(prev_hardL & prev_err),delta(prev_hardL & prev_err),0,[],0);

%% ------------------------------------------------------------------------
%% INDIVIDUAL MICE
%% ------------------------------------------------------------------------

%% overall alternation bias and by previous trial dificulty, for each animal
%  alternation is with respect to trial type (reward location)
%  also do n trials back

mice                    = unique(mouseid);
th.postRewBiasMice      = nan(1,numel(mice));
th.postEasyRewBiasMice  = nan(1,numel(mice));
th.postHardRewBiasMice  = nan(1,numel(mice));
th.postErrBiasMice      = nan(1,numel(mice));
th.postEasyRewBiasMice  = nan(1,numel(mice));
th.postHardRewBiasMice  = nan(1,numel(mice));
th.perfMice             = nan(1,numel(mice));
th.postRewPerfMice      = nan(1,numel(mice));
th.postErrPerfMice      = nan(1,numel(mice));

th.biasNbackMice        = zeros(numel(mice),cfg.nback);
th.biasNbackRewMice     = zeros(numel(mice),cfg.nback);
th.biasNbackErrMice     = zeros(numel(mice),cfg.nback);
th.biasConsecRewMice    = zeros(numel(mice),cfg.nback);
th.biasConsecErrMice    = zeros(numel(mice),cfg.nback);

for iMouse = 1:numel(mice)
  
  basel_rw  = psychometricFit(choice(mouseid == mice(iMouse) & prev_rw),delta(mouseid == mice(iMouse) & prev_rw),0,[],0);
  basel_er  = psychometricFit(choice(mouseid == mice(iMouse) & prev_err),delta(mouseid == mice(iMouse) & prev_err),0,[],0);
  
  rwL       = psychometricFit(choice((prev_easyL | prev_hardL) & prev_rw & mouseid == mice(iMouse)),delta((prev_easyL | prev_hardL) & prev_rw & mouseid == mice(iMouse)),0,[],0);
  rwR       = psychometricFit(choice((prev_easyR | prev_hardR) & prev_rw & mouseid == mice(iMouse)),delta((prev_easyR | prev_hardR) & prev_rw & mouseid == mice(iMouse)),0,[],0);
  easyRwL   = psychometricFit(choice(prev_easyL & prev_rw & mouseid == mice(iMouse)),delta(prev_easyL & prev_rw & mouseid == mice(iMouse)),0,[],0);
  hardRwL   = psychometricFit(choice(prev_hardL & prev_rw & mouseid == mice(iMouse)),delta(prev_hardL & prev_rw & mouseid == mice(iMouse)),0,[],0);
  easyRwR   = psychometricFit(choice(prev_easyR & prev_rw & mouseid == mice(iMouse)),delta(prev_easyR & prev_rw & mouseid == mice(iMouse)),0,[],0);
  hardRwR   = psychometricFit(choice(prev_hardR & prev_rw & mouseid == mice(iMouse)),delta(prev_hardR & prev_rw & mouseid == mice(iMouse)),0,[],0);
  
  errL      = psychometricFit(choice((prev_easyL | prev_hardL) & prev_err & mouseid == mice(iMouse)),delta((prev_easyL | prev_hardL) & prev_err & mouseid == mice(iMouse)),0,[],0);
  errR      = psychometricFit(choice((prev_easyR | prev_hardR) & prev_err & mouseid == mice(iMouse)),delta((prev_easyR | prev_hardR) & prev_err & mouseid == mice(iMouse)),0,[],0);
  easyErrL  = psychometricFit(choice(prev_easyL & prev_err & mouseid == mice(iMouse)),delta(prev_easyL & prev_err & mouseid == mice(iMouse)),0,[],0);
  hardErrL  = psychometricFit(choice(prev_hardL & prev_err & mouseid == mice(iMouse)),delta(prev_hardL & prev_err & mouseid == mice(iMouse)),0,[],0);
  easyErrR  = psychometricFit(choice(prev_easyR & prev_err & mouseid == mice(iMouse)),delta(prev_easyR & prev_err & mouseid == mice(iMouse)),0,[],0);
  hardErrR  = psychometricFit(choice(prev_hardR & prev_err & mouseid == mice(iMouse)),delta(prev_hardR & prev_err & mouseid == mice(iMouse)),0,[],0);
  
  th.postRewBiasMice(iMouse)      = 100.*(-nanmean(rwR.perfPsych - basel_rw.perfPsych)+nanmean(rwL.perfPsych - basel_rw.perfPsych));
  th.postEasyRewBiasMice(iMouse)  = 100.*(-nanmean(easyRwR.perfPsych - basel_rw.perfPsych)+nanmean(easyRwL.perfPsych - basel_rw.perfPsych));
  th.postHardRewBiasMice(iMouse)  = 100.*(-nanmean(hardRwR.perfPsych - basel_rw.perfPsych)+nanmean(hardRwL.perfPsych - basel_rw.perfPsych));
  th.postErrBiasMice(iMouse)      = -100.*(-nanmean(errR.perfPsych - basel_er.perfPsych)+nanmean(errL.perfPsych - basel_er.perfPsych));
  th.postEasyErrBiasMice(iMouse)  = -100.*(-nanmean(easyErrR.perfPsych - basel_er.perfPsych)+nanmean(easyErrL.perfPsych - basel_er.perfPsych));
  th.postHardErrBiasMice(iMouse)  = -100.*(-nanmean(hardErrR.perfPsych - basel_er.perfPsych)+nanmean(hardErrL.perfPsych - basel_er.perfPsych));
  
  th.perfMice(iMouse)             = 100.*sum(choice(mouseid == mice(iMouse)) == trialtype(mouseid == mice(iMouse)))./sum((mouseid == mice(iMouse)));
  th.postRewPerfMice(iMouse)      = 100.*sum(choice(mouseid == mice(iMouse) & prev_rw) == trialtype(mouseid == mice(iMouse) & prev_rw))./sum((prev_rw & mouseid == mice(iMouse)));
  th.postErrPerfMice(iMouse)      = 100.*sum(choice(mouseid == mice(iMouse) & prev_err) == trialtype(mouseid == mice(iMouse) & prev_err))./sum((prev_err & mouseid == mice(iMouse)));
  
  %% n trial back analysis
  deltaMouse  = delta(mouseid == mice(iMouse));
  rwMouse     = reward(mouseid == mice(iMouse));
  errMouse    = error(mouseid == mice(iMouse));
  chMouse     = choice(mouseid == mice(iMouse));
  
  for nBack = 1:cfg.nback
    
    delta_nB  = [zeros(1,nBack) deltaMouse(1:end-nBack)];
    rw_nB     = [zeros(1,nBack) rwMouse(1:end-nBack)];
    err_nB    = [zeros(1,nBack) errMouse(1:end-nBack)];
    
    rwL       = psychometricFit(chMouse(delta_nB<0 & rw_nB),deltaMouse(delta_nB<0 & rw_nB),0,[],0);
    rwR       = psychometricFit(chMouse(delta_nB>0 & rw_nB),deltaMouse(delta_nB>0 & rw_nB),0,[],0);
    errL      = psychometricFit(chMouse(delta_nB<0 & err_nB),deltaMouse(delta_nB<0 & err_nB),0,[],0);
    errR      = psychometricFit(chMouse(delta_nB>0 & err_nB),deltaMouse(delta_nB>0 & err_nB),0,[],0);
    
    th.biasNbackRewMice(iMouse,nBack) = 100.*(-nanmean(rwR.perfPsych - basel_rw.perfPsych) + nanmean(rwL.perfPsych - basel_rw.perfPsych));
    th.biasNbackErrMice(iMouse,nBack) = -100.*(-nanmean(errR.perfPsych- basel_er.perfPsych) + nanmean(errL.perfPsych- basel_er.perfPsych));
    
  end
  
  %% look at trial combinations, n consecutive rewards or errors on same side
  rw_i    = [0 rwMouse(1:end-1)];
  er_i    = [0 errMouse(1:end-1)];
  delta_i = sign([0 deltaMouse(1:end-1)]);
  delta_L = delta_i < 0;
  delta_R = delta_i > 0;
  
  for nBack = 1:cfg.nback
    
    % successive multiplications -- since vectors and 1's and 0's
    % multiplying 1's gives consecutive rewards or errors
    rw_i    = rw_i.*[zeros(1,nBack) rwMouse(1:end-nBack)];
    er_i    = er_i.*[zeros(1,nBack) errMouse(1:end-nBack)];
    delta_n = sign([zeros(1,nBack) deltaMouse(1:end-nBack)]);
    delta_i = delta_i .* delta_n; % compound delta sign vector (consec R or L)
    delta_j = delta_i > 0; % same sign multiplication is > 0
    
    rwL     = psychometricFit(chMouse(delta_L & delta_j & rw_i),deltaMouse(delta_L & delta_j & rw_i),0,[],0);
    rwR     = psychometricFit(chMouse(delta_R & delta_j & rw_i),deltaMouse(delta_R & delta_j & rw_i),0,[],0);
    errL    = psychometricFit(chMouse(delta_L & delta_j & er_i),deltaMouse(delta_L & delta_j & er_i),0,[],0);
    errR    = psychometricFit(chMouse(delta_R & delta_j & er_i),deltaMouse(delta_R & delta_j & er_i),0,[],0);
    
    th.biasConsecRewMice(iMouse,nBack) = 100.*(-nanmean(rwR.perfPsych-basel_rw.perfPsych)+nanmean(rwL.perfPsych-basel_rw.perfPsych));
    th.biasConsecErrMice(iMouse,nBack) = -100.*(-nanmean(errR.perfPsych-basel_er.perfPsych)+nanmean(errL.perfPsych-basel_er.perfPsych));
  
  end
end

%% sumarize mouse analysis
d_sem                      = sqrt(sum(~isnan(th.biasNbackRewMice(:,1)))-1);

th.postRewBiasEasyHard     = [nanmean(th.postEasyRewBiasMice) nanmean(th.postHardRewBiasMice)];
th.postRewBiasEasyHardSEM  = [nanstd(th.postEasyRewBiasMice) nanstd(th.postHardRewBiasMice)]./d_sem;
th.postErrBiasEasyHard     = [nanmean(th.postEasyErrBiasMice) nanmean(th.postHardErrBiasMice)];
th.postErrBiasEasyHardSEM  = [nanstd(th.postEasyErrBiasMice) nanstd(th.postHardErrBiasMice)]./d_sem;

th.biasNback               = nanmean(th.biasNbackMice);
th.biasNbackSEM            = nanstd(th.biasNbackMice)./d_sem;
th.biasNbackRew            = nanmean(th.biasNbackRewMice);
th.biasNbackRewSEM         = nanmean(th.biasNbackRewMice)./d_sem;
th.biasNbackErr            = nanmean(th.biasNbackErrMice);
th.biasNbackErrSEM         = nanmean(th.biasNbackErrMice)./d_sem;

th.biasConsecRew           = nanmean(th.biasConsecRewMice);
th.biasConsecRewSEM        = nanmean(th.biasConsecRewMice)./d_sem;
th.biasConsecErr           = nanmean(th.biasConsecErrMice);
th.biasConsecErrSEM        = nanmean(th.biasConsecErrMice)./d_sem;

th.stats                   = getStats(th);

%% calculate psychometrics removing perseverence blocks
goodTrials        = detectPerseverenceBouts(lg,cfg.nConsecPers);
th.psychNoBouts   = psychometricFit(choice(goodTrials),delta(goodTrials),0,[],0);

end

%% get summary statistics
function stats = getStats(th)

% correlation between post-prev_rw and post-prev_err bias for individual mice
validx = find(~isnan(th.postRewBiasMice) & ~isnan(th.postErrBiasMice));
[stats.corr_postRW_vs_postErr,stats.p_postRW_vs_postErr] = ...
  corr(abs(th.postRewBiasMice(validx))',abs(th.postErrBiasMice(validx))');

% significance of post-prev_rw and post-prev_err bias at population level
if lillietest(th.postRewBiasMice)
  stats.postRewBias_pval     = signrank(th.postRewBiasMice);
  stats.postRewBias_test     = 'signrank';
else
  [~,stats.postRewBias_pval] = ttest(th.postRewBiasMice);
  stats.postRewBias_test     = 'ttest';
end
if lillietest(th.postErrBiasMice)
  stats.postErrBias_pval     = signrank(th.postErrBiasMice);
  stats.postErrBias_test     = 'signrank';
else
  [~,stats.postErrBias_pval] = ttest(th.postErrBiasMice);
  stats.postErrBias_test     = 'ttest';
end

% post-prev_rw ~= post-prev_err?
if lillietest(th.postRewBiasMice) || lillietest(th.postErrBiasMice)
  stats.postRewVSpostErrBias_pval     = signrank(abs(th.postRewBiasMice),abs(th.postErrBiasMice));
  stats.postRewVSpostErrBias_test     = 'signrank';
else
  [~,stats.postRewVSpostErrBias_pval] = ttest(abs(th.postRewBiasMice),abs(th.postErrBiasMice));
  stats.postRewVSpostErrBias_test     = 'ttest';
end

% effect of previous trial difficulty: 2-way RM anova
nmice   = numel(th.postEasyRewBiasMice);
data    = [th.postEasyRewBiasMice'; ...
           th.postHardRewBiasMice'; ...
           th.postEasyErrBiasMice'; ...
           th.postHardErrBiasMice'];
mvec    = repmat((1:nmice)',[4 1]);
outvec  = [ones(nmice*2,1); -ones(nmice*2,1)];
diffvec = repmat([ones(nmice,1); 2.*ones(nmice,1)],[2 1]);

stats.trialDifficulty_testname = '2-way RM ANOVA';
stats.trialDifficulty_varname  = {'outcome','difficulty','mouseid'};
[stats.trialDifficulty_pvals,stats.trialDifficulty_table]   =       ...
                anovan(data,{outvec, diffvec, mvec},                ...
                      'varnames',{'outcome','difficulty','mouseid'},...
                      'display','off');

end
