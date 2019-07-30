function rc = logRegressionFromConcatLog(lg,minNT,nBins)

% rc = logRegressionFromConcatLog(lg,minNT)
% quantifies per-mouse logistic regression. lg is flattened log, minNT is
% minimum number of trials to be included in the analysis

if nargin < 2
  minNT    = 1000;
end
if nargin < 3
  nBins    = 5;
end

lg         = cleanupConcatLog(lg,minNT);
mice       = unique(lg.mouseID);
nmice      = numel(mice);

rc.overall = revCorr_logistic(lg.trialType,lg.choice,lg.cuePos_R,lg.cuePos_L,[],0,nBins);
rc.bins    = rc.overall.bins;
rc.values  = rc.overall.values;

% stats per mouse
rc.goodMiceID            = analysisParams.miceBehav(mice);
rc.ntrials               = zeros(1,nmice);
rc.goodValues            = zeros(numel(rc.values),nmice);
rc.slopes                = zeros(1,nmice);
rc.intercepts            = zeros(1,nmice);
rc.slopesPoly1           = zeros(1,nmice);
rc.interceptsPoly1       = zeros(1,nmice);
rc.quadratic             = zeros(1,nmice);
rc.peaks                 = zeros(1,nmice);
rc.isSigSlope            = zeros(1,nmice);
rc.isSigSlopeBoot        = zeros(1,nmice);
rc.isSigQuadr            = zeros(1,nmice);
rc.nSessGoodMice         = zeros(1,nmice);
rc.r2lin                 = zeros(1,nmice);
rc.r2exp                 = zeros(1,nmice);
rc.expFactor             = zeros(1,nmice);
rc.derivQuadr            = zeros(1,nmice);
rc.decayIndex            = zeros(1,nmice);
rc.isSigDerivQuadrBoot   = zeros(1,nmice);
rc.isSigDecayIndexBoot   = zeros(1,nmice);


% for individual mice
for mm = 1:numel(mice)
  [choice, cpr, cpl, tt]    ...
                 = selectMouseTrials(lg, mice(mm), 'choice', 'cuePos_R', 'cuePos_L', 'trialType');
  rc.mice(mm)    = revCorr_logistic(tt, choice, cpr, cpl ,[],0,nBins,true); 
  rc.ntrials(mm) = numel(choice);
  
  % general info
  rc.nSessGoodMice(:,mm) = numel(unique(lg.sessionID(lg.mouseID == mice(mm))));
  x                      = rc.bins(1:end-1)'+mode(diff(rc.bins))/2;
  y                      = rc.mice(mm).values;
  rc.decayIndex(:,mm)    = rc.mice(mm).decayIndex;%nanmean(y(end-1:end))/nanmean(y(1:2));
  rc.isSigDecayIndexBoot(:,mm) = rc.mice(mm).decayIndex_p<.05;
  rc.goodValues(:,mm)    = rc.mice(mm).values;
  rc.peaks(:,mm)         = x(y==max(y));
  
  % quadratic fit
  [temp,stats]            = fit(x,y,'poly2');
  yhat                    = fiteval(temp,x);
  rc.derivQuadr(:,mm)     = nanmean(diff(yhat));
  coeffs                  = coeffvalues(temp);
  ci                      = confint(temp);
  rc.r2poly2(mm)          = stats.rsquare;
  rc.fitStats_poly2(mm)   = stats;
  rc.slopes(:,mm)         = coeffs(2);
  rc.intercepts(:,mm)     = coeffs(3);
  rc.quadratic(:,mm)      = coeffs(1);
  if ci(1,2) <= 0 && ci(2,2) >= 0; sig = 0; else sig = 1; end
  if ci(1,1) <= 0 && ci(2,1) >= 0; sig2 = 0; else sig2 = 1; end
  rc.isSigSlope(:,mm)     = sig;
  rc.isSigQuadr(:,mm)     = sig2;
  
  % linear fit
  [temp,stats]                = fit(x,y,'poly1');
  coeffs                      = coeffvalues(temp);
  rc.r2lin(mm)                = stats.rsquare;
  rc.fitStats_lin(mm)         = stats;
  rc.slopesPoly1(:,mm)        = coeffs(1);
  rc.interceptsPoly1(:,mm)    = coeffs(2);
  
  % exponential fit
  [temp,stats]                = fit(rc.bins(1:end-1)'+mode(diff(rc.bins))/2,rc.mice(mm).values,...
    fittype(@(a,b,c,x) a+b*exp(c*x)),'startPoint',[.1 0 .2],'lower',[0 -5 -1],'upper',[1 5 1]);
  temp                        = coeffvalues(temp);
  rc.r2exp(mm)                = stats.rsquare;
  rc.fitStats_exp(mm)         = stats;
  rc.expFactor(:,mm)          = temp(3);
end

rc.nGoodMice  = size(rc.goodValues,2);
rc.sem        = nanstd(rc.goodValues,0,2)./sqrt(rc.nGoodMice-1);

for mm = 1:numel(rc.mice)
  rc.mice(mm).bootStd = rc.mice(mm).sem;
end