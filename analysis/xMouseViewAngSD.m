function vsd = xMouseViewAngSD(lg,minNumTrials,posBins)

% vsd = xMouseViewAngSD(lg,minNumTrials,posBins)
% average per-mouse view angle trajectory standard deviation
% lg is flattened behavior log
% minNumTrials is minimal trial count to be included in analysis (default 0)
% posBins is vector defining binning of y positions in cm, default 0:5:300


if nargin < 2; minNumTrials   = 0;          end
if nargin < 3; posBins        = 0:5:300;    end

lg    = cleanupConcatLog(lg, minNumTrials);

vsd.mice         = unique(lg.mouseID);
vsd.nmice        = length(vsd.mice);
vsd.viewAngSD    = zeros(1,vsd.nmice);
vsd.goodSessPerc = nan(1,vsd.nmice);

for iMouse = 1:vsd.nmice
  pos     = selectMouseTrials(lg, vsd.mice(iMouse), 'pos');
  viewAng = sampleViewAngleVsY(pos, posBins, false);

  % average std of view angle trajectories  
  vsd.viewAngSD(iMouse)= mean(std(viewAng));
end

vsd.mouseAvg         = mean(vsd.viewAngSD);
vsd.mouseSem         = std(vsd.viewAngSD)./sqrt(vsd.nmice-1);
