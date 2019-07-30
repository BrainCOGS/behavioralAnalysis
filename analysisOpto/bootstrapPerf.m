function boot = bootstrapPerf(choice,trialType,laserON,nCues_RminusL,mouseID,sessionID,niter,drawMiceSessions,doFit,trialDiff)

% boot = bootstrapPerf(choice,trialType,laserON,nCues_RminusL,mouseID,sessionID,niter,drawMiceSessions,doFit)
%
% svoboda randomly samples session, mouse and trials. implement this later:
% We tested against the null hypothesis that each photoinhibition site did not cause a performance change.
% Performance changes must be interpreted against behavioral variability. We performed bootstrap to consider
% the variability across mice, sessions, and trials (Efron and Tibshirani, 1994). For each cycle of bootstrap,
% repeated 106 times, we randomly sampled with replacement (1) animals, (2) sessions performed by each animal,
% and (3) trials within each session. We computed control performance and performance under photoinhibition
% for each condition (i.e., ??lick right?? and ??lick left?? trial per- formances for each cortical location).
% The p value for each photoinhibition condition was the fraction of times the performance change from
% the control condition changed sign (if photoinhibition showed a mean decrease in perfor- mance from the control,
% p value for this condition was the number of times it showed an increase in performance during bootstrap).
% drawMiceSessions is a flag that if set will repeat the procedure
% described above, otherwise bootstrapping iterations draw blindly to mice
% and sessions
%
% OUTPUT:
%   boot.ctrl.sd_percCorrect
%   boot.lsr.sd_percCorrect
%   boot.ctrl.sd_percFinish
%   boot.lsr.sd_percFinish
%   boot.ctrl.sd_biasPsych
%   boot.lsr.sd_biasPsych
%
%   boot.p_percCorrect
%   boot.p_percFinish
%   boot.p_biasPsych
%
% lp feb 2016

%% defualts
rng('default')

if nargin < 9;  niter            = 5000; end
if nargin < 10; drawMiceSessions = true; end
if nargin < 11; doFit            = true; end
if nargin < 12; trialDiff        = []; end

%% calculate overall performance for comparison at each bootstrap iteration
pcCtrl  = 100*sum(trialType(laserON == 0) == choice(laserON == 0))./sum(laserON==0);
pcLsr   = 100*sum(trialType(laserON == 1) == choice(laserON == 1))./sum(laserON==1);
pcLsrNm = (pcLsr - pcCtrl) / (pcCtrl - 50);
pcCtrlR = 100*sum(trialType(laserON == 0 & trialType == analysisParams.rightCode) == choice(laserON == 0 & trialType == analysisParams.rightCode))./sum(laserON==0 & trialType == analysisParams.rightCode);
pcLsrR  = 100*sum(trialType(laserON == 1 & trialType == analysisParams.rightCode) == choice(laserON == 1 & trialType == analysisParams.rightCode))./sum(laserON==1 & trialType == analysisParams.rightCode);
pcCtrlL = 100*sum(trialType(laserON == 0 & trialType == analysisParams.leftCode) == choice(laserON == 0 & trialType == analysisParams.leftCode))./sum(laserON==0 & trialType == analysisParams.leftCode);
pcLsrL  = 100*sum(trialType(laserON == 1 & trialType == analysisParams.leftCode) == choice(laserON == 1 & trialType == analysisParams.leftCode))./sum(laserON==1 & trialType == analysisParams.leftCode);
pfCtrl  = 100*sum(choice(laserON == 0)~=analysisParams.nilCode)./sum(laserON==0);
pfLsr   = 100*sum(choice(laserON == 1)~=analysisParams.nilCode)./sum(laserON==1);
bctrl   = pcCtrlR - pcCtrlL;
blsr    = pcLsrR  - pcLsrL;
bactrl  = abs(pcCtrlR - pcCtrlL);
balsr   = abs(pcLsrR  - pcLsrL);
lpsCtrl = 100*sum(trialType(laserON == 0 & abs(nCues_RminusL) >= 10) ~= choice(laserON == 0 & abs(nCues_RminusL) >= 10))./sum(laserON==0 & abs(nCues_RminusL) >= 10);
lpsLsr  = 100*sum(trialType(laserON == 1 & abs(nCues_RminusL) >= 10) ~= choice(laserON == 1 & abs(nCues_RminusL) >= 10))./sum(laserON==1 & abs(nCues_RminusL) >= 10);

if isempty(trialDiff); trialDiff = -abs(nCues_RminusL); end
iseasy       = trialDiff < median(trialDiff);
ishard       = trialDiff > median(trialDiff);
pcCtrl_easy  = 100*sum(trialType(laserON == 0 & iseasy) == choice(laserON == 0 & iseasy))./sum(laserON==0 & iseasy);
pcLsr_easy   = 100*sum(trialType(laserON == 1 & iseasy) == choice(laserON == 1 & iseasy))./sum(laserON==1 & iseasy);
pcCtrl_hard  = 100*sum(trialType(laserON == 0 & ishard) == choice(laserON == 0 & ishard))./sum(laserON==0 & ishard);
pcLsr_hard   = 100*sum(trialType(laserON == 1 & ishard) == choice(laserON == 1 & ishard))./sum(laserON==1 & ishard);

%% for psychometric curve (biasPsych % slope estimation)
if doFit
  psychCtrl   = psychometricFit(choice(laserON == 0),nCues_RminusL(laserON == 0),false,[],true);
  psychLsr    = psychometricFit(choice(laserON == 1),nCues_RminusL(laserON == 1),false,[],true);
  sbpsych     = sign(mean(psychLsr.perfPsychJ_binned-psychCtrl.perfPsychJ_binned));
end

% sign of differences between lsr and ctrl for significance calculation
spc    = sign(pcLsr-pcCtrl);
spcNm  = sign(pcLsrNm);
spcR   = sign(pcLsrR-pcCtrlR);
spcL   = sign(pcLsrL-pcCtrlL);
spcE   = sign(pcLsr_easy-pcCtrl_easy);
spcH   = sign(pcLsr_hard-pcCtrl_hard);
spf    = sign(pfLsr-pfCtrl);
sb     = sign(blsr-bctrl);
sba    = sign(balsr-bactrl);
slps   = sign(lpsLsr - lpsCtrl);

if doFit
  if ~isempty(psychLsr.fit) && ~isempty(psychCtrl.fit)
    slp     = sign(psychLsr.fit.slope - psychCtrl.fit.slope);
    fitFlag = 1;
  else
    slp     = nan;
    fitFlag = 0;
  end
else
  slp       = nan;
end

%% bootstrap
pcCtrl    = zeros(1,niter);
pcCtrlR   = zeros(1,niter);
pcCtrlL   = zeros(1,niter);
pcCtrlE   = zeros(1,niter);
pcCtrlH   = zeros(1,niter);
pfCtrl    = zeros(1,niter);
bctrl     = zeros(1,niter);
bactrl    = zeros(1,niter);
slpCtrl   = zeros(1,niter);
pcLsr     = zeros(1,niter);
pcLsrNm   = zeros(1,niter);
pcLsrR    = zeros(1,niter);
pcLsrL    = zeros(1,niter);
pcLsrE    = zeros(1,niter);
pcLsrH    = zeros(1,niter);
pfLsr     = zeros(1,niter);
blsr      = zeros(1,niter);
balsr     = zeros(1,niter);
slpLsr    = zeros(1,niter);
biasPsych = zeros(1,niter);
ntrials   = length(choice);
lpsCtrl   = zeros(1,niter);
lpsLsr    = zeros(1,niter);

for ii = 1:niter
  if drawMiceSessions % bootstrap over mice, sessions & trials
    mlist       = unique(mouseID);
    nmice       = length(mlist);
    midx        = mlist(randsample(nmice,nmice,true)); % sample with replacement from mice
    
    % have to loop here to ensure unique sessions for each mouse
    idx = [];
    for mm = 1:nmice
      tidx    = find(mouseID == midx(mm));
      nsess   = length(unique(sessionID(tidx)));
      sidx    = randsample(nsess,nsess,true); % sample with replacement from sessions
      for ss = 1:nsess
        thisl = sum(sessionID(tidx)==sidx(ss));
        idx(end+1:end+thisl) = tidx(sessionID(tidx)==sidx(ss));
      end
    end
    idx         = idx(randsample(length(idx),length(idx),true)); % sample with replacement from trials
  else                % or bootsrap over concatenated trials regardless of mouse and session
    idx         = randsample(ntrials,ntrials,true); % sample with replacement from all trials
  end
  
  tchoice         = choice(idx);
  ttrialType      = trialType(idx);
  tlaserON        = laserON(idx);
  tnCues_RminusL  = nCues_RminusL(idx);
  tiseasy         = iseasy(idx);
  tishard         = ishard(idx);
  
  % on rare occasions some mice don't have a given location
  if sum(tlaserON) == 0
    pcCtrl(ii)    = nan;
    pcLsr(ii)     = nan;
    pcLsrNm(ii)   = nan;
    pcCtrlR(ii)   = nan;
    pcLsrR(ii)    = nan;
    pcCtrlL(ii)   = nan;
    pcLsrL(ii)    = nan;
    pcCtrlE(ii)   = nan;
    pcLsrE(ii)    = nan;
    pcCtrlH(ii)   = nan;
    pcLsrH(ii)    = nan;
    pfCtrl(ii)    = nan;
    pfLsr(ii)     = nan;
    bctrl(ii)     = nan;
    blsr(ii)      = nan;
    bactrl(ii)    = nan;
    balsr(ii)     = nan;
    biasPsych(ii) = nan;
    slpCtrl(ii)   = nan;
    slpLsr(ii)    = nan;
    lpsCtrl(ii)   = nan;
    lpsLsr(ii)    = nan;
  else
    
    % calculate indicators for this iteration
    pcCtrl(ii)  = 100*sum(ttrialType(tlaserON == 0) == tchoice(tlaserON == 0))./sum(tlaserON==0);
    pcLsr(ii)   = 100*sum(ttrialType(tlaserON == 1) == tchoice(tlaserON == 1))./sum(tlaserON==1);
    pcLsrNm(ii) = (pcLsr(ii) - pcCtrl(ii)) / (pcCtrl(ii) - 50);
    pcCtrlE(ii) = 100*sum(ttrialType(tlaserON == 0 & tiseasy) == tchoice(tlaserON == 0 & tiseasy))./sum(tlaserON==0 & tiseasy);
    pcLsrE(ii)  = 100*sum(ttrialType(tlaserON == 1 & tiseasy) == tchoice(tlaserON == 1 & tiseasy))./sum(tlaserON==1 & tiseasy);
    pcCtrlH(ii) = 100*sum(ttrialType(tlaserON == 0 & tishard) == tchoice(tlaserON == 0 & tishard))./sum(tlaserON==0 & tishard);
    pcLsrH(ii)  = 100*sum(ttrialType(tlaserON == 1 & tishard) == tchoice(tlaserON == 1 & tishard))./sum(tlaserON==1 & tishard);
    pcCtrlR(ii) = 100*sum(ttrialType(tlaserON == 0 & ttrialType == analysisParams.rightCode) == tchoice(tlaserON == 0 & ttrialType == analysisParams.rightCode))./sum(tlaserON==0 & ttrialType == analysisParams.rightCode);
    pcLsrR(ii)  = 100*sum(ttrialType(tlaserON == 1 & ttrialType == analysisParams.rightCode) == tchoice(tlaserON == 1 & ttrialType == analysisParams.rightCode))./sum(tlaserON==1 & ttrialType == analysisParams.rightCode);
    pcCtrlL(ii) = 100*sum(ttrialType(tlaserON == 0 & ttrialType == analysisParams.leftCode) == tchoice(tlaserON == 0 & ttrialType == analysisParams.leftCode))./sum(tlaserON==0 & ttrialType == analysisParams.leftCode);
    pcLsrL(ii)  = 100*sum(ttrialType(tlaserON == 1 & ttrialType == analysisParams.leftCode) == tchoice(tlaserON == 1 & ttrialType == analysisParams.leftCode))./sum(tlaserON==1 & ttrialType == analysisParams.leftCode);
    pfCtrl(ii)  = 100*sum(tchoice(tlaserON == 0)~=analysisParams.nilCode)./sum(tlaserON==0);
    pfLsr(ii)   = 100*sum(tchoice(tlaserON == 1)~=analysisParams.nilCode)./sum(tlaserON==1);
    bctrl(ii)   = pcCtrlR(ii)-pcCtrlL(ii);
    blsr(ii)    = pcLsrR(ii)-pcLsrL(ii);
    bactrl(ii)  = abs(pcCtrlR(ii)-pcCtrlL(ii));
    balsr(ii)   = abs(pcLsrR(ii)-pcLsrL(ii));
    lpsCtrl(ii) = 100*sum(ttrialType(tlaserON == 0 & abs(tnCues_RminusL) >= 10) ~= tchoice(tlaserON == 0 & abs(tnCues_RminusL) >= 10))./sum(tlaserON==0 & abs(tnCues_RminusL) >= 10);
    lpsLsr(ii)  = 100*sum(ttrialType(tlaserON == 1 & abs(tnCues_RminusL) >= 10) ~= tchoice(tlaserON == 1 & abs(tnCues_RminusL) >= 10))./sum(tlaserON==1 & abs(tnCues_RminusL) >= 10);

    
    if doFit
      if fitFlag
        psychCtrl       = psychometricFit(tchoice(tlaserON == 0),tnCues_RminusL(tlaserON == 0),false,[],true);
        psychLsr        = psychometricFit(tchoice(tlaserON == 1),tnCues_RminusL(tlaserON == 1),false,[],true);
        biasPsych(ii)   = mean(psychLsr.perfPsychJ_binned-psychCtrl.perfPsychJ_binned);
        try
          slpCtrl(ii)   = psychCtrl.fit.slope;
        catch
          slpCtrl(ii)   = nan;
        end
        try
          slpLsr(ii)    = psychLsr.fit.slope;
        catch
          slpLsr(ii)    = nan;
        end
      else
        biasPsych(ii)   = nan;
        slpCtrl(ii)     = nan;
        slpLsr(ii)      = nan;
      end
    else
      biasPsych(ii)    = nan;
      slpCtrl(ii)      = nan;
      slpLsr(ii)       = nan;
    end
  end
end

% calculate p-value, see svoboda explanation above
boot.ctrl.sd_percCorrect        = nanstd(pcCtrl);
boot.lsr.sd_percCorrect         = nanstd(pcLsr);
boot.lsr.sd_percCorrect_norm    = nanstd(pcLsrNm);
boot.ctrl.sd_percCorrect_easy   = nanstd(pcCtrlE);
boot.lsr.sd_percCorrect_easy    = nanstd(pcLsrE);
boot.ctrl.sd_percCorrect_hard   = nanstd(pcCtrlH);
boot.lsr.sd_percCorrect_hard    = nanstd(pcLsrH);
boot.ctrl.sd_percCorrect_R      = nanstd(pcCtrlR);
boot.lsr.sd_percCorrect_R       = nanstd(pcLsrR);
boot.ctrl.sd_percCorrect_L      = nanstd(pcCtrlL);
boot.lsr.sd_percCorrect_L       = nanstd(pcLsrL);
boot.ctrl.sd_percFinish         = nanstd(pfCtrl);
boot.lsr.sd_percFinish          = nanstd(pfLsr);
boot.ctrl.sd_bias               = nanstd(bctrl);
boot.lsr.sd_bias                = nanstd(blsr);
boot.ctrl.sd_bias_abs           = nanstd(bactrl);
boot.lsr.sd_bias_abs            = nanstd(balsr);
boot.ctrl.sd_slope              = nanstd(slpCtrl);
boot.lsr.sd_slope               = nanstd(slpLsr);
boot.ctrl.sd_lapse              = nanstd(lpsCtrl);
boot.lsr.sd_lapse               = nanstd(lpsLsr);

boot.lsr.diff_percCorrect_norm  = nanmean(pcLsrNm);
boot.diff_percCorrect_mean      = nanmean(pcLsr - pcCtrl);
boot.diff_percCorrect_std       = nanstd(pcLsr - pcCtrl);
boot.diff_percCorrect_easy_mean = nanmean(pcLsrE - pcCtrlE);
boot.diff_percCorrect_easy_std  = nanstd(pcLsrE - pcCtrlE);
boot.diff_percCorrect_hard_mean = nanmean(pcLsrH - pcCtrlH);
boot.diff_percCorrect_hard_std  = nanstd(pcLsrH - pcCtrlH);
boot.diff_percCorrect_R_mean    = nanmean(pcLsrR - pcCtrlR);
boot.diff_percCorrect_R_std     = nanstd(pcLsrR - pcCtrlR);
boot.diff_percCorrect_L_mean    = nanmean(pcLsrL - pcCtrlL);
boot.diff_percCorrect_L_std     = nanstd(pcLsrL - pcCtrlL);
boot.diff_percFinish_mean       = nanmean(pfLsr - pfCtrl);
boot.diff_percFinish_std        = nanstd(pfLsr - pfCtrl);
boot.diff_bias_mean             = nanmean(blsr - bctrl);
boot.diff_bias_std              = nanstd(blsr - bctrl);
boot.diff_bias_abs_mean         = nanmean(balsr - bactrl);
boot.diff_bias_abs_std          = nanstd(balsr - bactrl);
boot.diff_slope_mean            = nanmean(slpLsr - slpCtrl);
boot.diff_slope_std             = nanstd(slpLsr - slpCtrl);
boot.diff_lapse_mean            = nanmean(lpsLsr - lpsCtrl);
boot.diff_lapse_std             = nanstd(lpsLsr - lpsCtrl);

boot.p_percCorrect              = nansum(sign(pcLsr-pcCtrl)~=spc)/niter;
boot.p_percCorrect_norm         = nansum(sign(pcLsrNm)~=spcNm)/niter;
boot.p_percCorrect_easy         = nansum(sign(pcLsrE-pcCtrlE)~=spcE)/niter;
boot.p_percCorrect_hard         = nansum(sign(pcLsrH-pcCtrlH)~=spcH)/niter;
boot.p_percCorrect_R            = nansum(sign(pcLsrR-pcCtrlR)~=spcR)/niter;
boot.p_percCorrect_L            = nansum(sign(pcLsrL-pcCtrlL)~=spcL)/niter;
boot.p_percFinish               = nansum(sign(pfLsr-pfCtrl)~=spf)/niter;
boot.p_bias                     = nansum(sign(blsr-bctrl)~=sb)/niter;
boot.p_bias_abs                 = nansum(sign(balsr-bactrl)~=sba)/niter;
boot.p_slope                    = nansum(sign(slpLsr-slpCtrl)~=slp)/niter;
boot.p_lapse                    = nansum(sign(lpsLsr-lpsCtrl)~=slps)/niter;

boot.biasPsych                  = nanmean(biasPsych);
boot.sd_biasPsych               = nanstd(biasPsych);
boot.p_biasPsych                = nansum(sign(biasPsych)~=sbpsych)/niter;

