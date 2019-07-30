function motorErrs = xMouseMotorErrors(lg,whichType)

% motorErrs = xMouseMotorErrors(lg,whichType)
% average per-mouse motor errors 
% lg is flattened behavior log
% which type is string speficying the type of motor event, default 'any'
% refer to unusualMotorEvents() for a full list of options

if nargin < 2; whichType = 'any'; end

mice        = unique(lg.mouseID);
motorErrs   = nan(numel(mice), 1);

for iMouse = 1:numel(mice)
  trialID            = unusualMotorEvents(lg,mice(iMouse),whichType);
  motorErrs(iMouse)  = sum(trialID)/numel(trialID) * 100;
end