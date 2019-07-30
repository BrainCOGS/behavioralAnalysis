function [vAngDecodeAcc,vAng] = xMouseViewAngleDecode(lg,posBins)

% [vAngDecodeAcc,vAng] = xMouseViewAngleDecode(lg,posBins)
% average per-mouse choice decoding based on view angles
% lg is flattened behavior log
% posBins is vector defining binning of y positions in cm, default 0:300

if nargin < 2
  posBins = 0:300;
end

mice             = unique(lg.mouseID);
vAngDecodeAcc    = nan(numel(posBins), numel(mice));
viewAngle        = sampleViewAngleVsY(lg.pos, posBins);

for iMouse = 1:numel(mice)
  [choice,sel]          = selectMouseTrials(lg, mice(iMouse), 'choice');
  viewAng               = viewAngle(:,sel);
  viewAng               = viewAng + (rand(size(viewAng)) - 0.5)*1e-4;     % prevent decoding degeneracies
  choiceVAng(:,iMouse)  = arrayfun(@(x) viewAng(:,choice == x), 0:1, 'UniformOutput', false);
  for iPos = 1:numel(posBins)
    try
      [threshold, decodeFrac]     = equalCDFThreshold(choiceVAng{1,iMouse}(iPos,:), choiceVAng{2,iMouse}(iPos,:), [], [], false);
      vAngDecodeAcc(iPos,iMouse)  = mean(decodeFrac);
    catch
      vAngDecodeAcc(iPos,iMouse)  = nan;
    end
  end
end

vAng.decodeAcc     = vAngDecodeAcc;
vAng.posBins       = posBins;
vAng.decodeAcc_avg = nanmean(vAngDecodeAcc,2);
vAng.decodeAcc_sem = nanstd(vAngDecodeAcc,0,2)./sqrt(numel(mice)-1);