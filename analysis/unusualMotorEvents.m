function [trialID, delta, choice] = unusualMotorEvents(lg,mouseID,eventType)

% [trialID,delta] = unusualMotorEvents(lg,mouseID,eventType)
%
% INPUT: lg is flattened structure, mouseID (number, def. -1 metamouse),
%  eventType is string:
%      'cueTurnAround': view angles sharper than 60deg during cue period
%       'excessTravel': regular 10% maze length excess threshold
%           'lowSpeed': calculated for individual animals, bottom 10% of trial
%          'earlyTurn': unusually large angles between 250 and 295 cm, with sign changes
% 'enteredOppositeArm': went into one arm and ended up in the other
%
% OUTPUT:
% trialID is logical vector for trials with unual motor events:
% delta is just #R - #L, choice as usual, for convenience
% 
% Lucas Pinto (lpinto@princeton.edu)

if nargin < 2 || isempty(mouseID)
  mouseID   = -1; % metamouse
end
if nargin < 3
  eventType = 'any';
end

[pos, speed, excessTravel, delta, choice, mid] ...
                    = selectMouseTrials(lg, mouseID, 'pos', 'speedStem', 'excessTravel', 'nCues_RminusL', 'choice', 'mouseID');
xpos                = cellfun(@(x) x(:,1), pos, 'UniformOutput', false);
ypos                = cellfun(@(x) x(:,2), pos, 'UniformOutput', false);
vang                = cellfun(@(x) x(:,3), pos, 'UniformOutput', false);

switch eventType
  case 'cueTurnAround' % view angles sharper than 80deg during cue
    trialID    = cellfun(@(x,y) any(x<200 & abs(y)>60), ypos, vang, 'UniformOutput', false);
    trialID    = [trialID{:}];
    
  case 'excessTravel' % regular 10% maze length excess threshold
    trialID    = excessTravel > .1;
    
  case 'lowSpeed' % calculated for individual animals, bottom 10% of trials
    if mouseID ~= -1
      trialID  = speed < prctile(speed,10);
    else
      trialID  = false(size(speed));
      mice     = unique(mid);
      for iMouse = 1:numel(mice)
        [ispeed, isel] = selectMouseTrials(lg, mice(iMouse), 'speedStem');
        isel           = find(isel == 1);
        trialID(isel(ispeed < prctile(ispeed,10))) = true;
      end
      
    end
    
  case 'earlyTurn' % unusually large angles between 250 and 295 cm, with sign changes
    va              = cellfun(@(x,y) x(y > 250 & y < 295), vang, ypos, 'UniformOutput', false);
    badva           = cellfun(@isempty, va, 'UniformOutput', false);
    badva           = [badva{:}];
    trialID         = false(size(badva));
    trials          = cellfun(@(x) any(abs(x)  > abs(x(1))+30)                    & ...
                              any(abs(x(find(abs(x)>abs(x(1))+30,1,'first'):end)) < ...
                                  abs(x(find(abs(x)>abs(x(1))+30,1,'first')))-2),   ...
                              va(~badva), 'UniformOutput', false);
    trialID(~badva) = [trials{:}];
    
  case 'enteredOppositeArm' % went into one arm and ended up in the other
    xp              = cellfun(@(x,y) x(y > 300), xpos, ypos, 'UniformOutput', false);
    badxp           = cellfun(@isempty, xp, 'UniformOutput', false);
    badxp           = [badxp{:}];
    trialID         = false(size(badxp));
    trials          = cellfun(@(x) (x(end) > 0 & any(x < -5)) | (x(end) < 0 & any(x > 5)), ...
                              xp(~badxp), 'UniformOutput', false);
    trialID(~badxp) = [trials{:}];
    
  case 'any'
    
    trialsc         = cellfun(@(x,y) any(x<200 & abs(y)>60), ypos, vang, 'UniformOutput', false);
    trialsc         = [trialsc{:}];
    
    va              = cellfun(@(x,y) x(y > 250 & y < 295), vang, ypos, 'UniformOutput', false);
    badva           = cellfun(@isempty, va, 'UniformOutput', false);
    badva           = [badva{:}];
    trialsi         = false(size(badva));
    trials          = cellfun(@(x) any(abs(x)  > abs(x(1))+30)                    & ...
                              any(abs(x(find(abs(x)>abs(x(1))+30,1,'first'):end)) < ...
                                  abs(x(find(abs(x)>abs(x(1))+30,1,'first')))-2),   ...
                              va(~badva), 'UniformOutput', false);
    trialsi(~badva) = [trials{:}];
    
    xp              = cellfun(@(x,y) x(y > 300), xpos, ypos, 'UniformOutput', false);
    badxp           = cellfun(@isempty, xp, 'UniformOutput', false);
    badxp           = [badxp{:}];
    trialsx         = false(size(badxp));
    trials          = cellfun(@(x) (x(end) > 0 & any(x < -5)) | (x(end) < 0 & any(x > 5)), ...
                              xp(~badxp), 'UniformOutput', false);
    trialsx(~badxp) = [trials{:}];
    
    trialID         = excessTravel > .1          | ...
                      speed < prctile(speed,10)  | ...
                      trialsi                    | ...
                      trialsx                    | ...
                      trialsc;
                            
    
end
