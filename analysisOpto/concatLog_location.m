%% flatten log structure and select relevant trials
function log = concatLog_location(logs,filters,lidx,perfTh)

vectorNames = {'laserON','galvoPosIdxMaster','trialType','choice','currMaze', ...
               'nCues_RminusL','excessTravel','speedStem','trialDur'        , ...
               'cuePos_R','cuePos_L','pos','meanPerfBlockCtrl','stimID'     , ...
               'nCues_R','nCues_L'};
           
for ii = 1:length(vectorNames)
    eval(sprintf('log.%s = [logs(:).%s];',vectorNames{ii},vectorNames{ii}));
end

% create mouse and session vectors for bootstrapping
log.mouseID = []; log.sessionID = []; sessct = zeros(1,length(analysisParams.mice));
for ii = 1:length(logs)
    mid = find(strcmp(analysisParams.mice,logs(ii).info.mouseID));
    log.mouseID(end+1:end+length(logs(ii).trialType)) = ...
        ones(1,length(logs(ii).trialType)).*mid;
    sessct(mid) = sessct(mid) + 1;
    log.sessionID(end+1:end+length(logs(ii).trialType)) = ...
        ones(1,length(logs(ii).trialType)).*sessct(mid);
end

vectorNames{end+1} = 'mouseID';
vectorNames{end+1} = 'sessionID';

% threshold performance?
if ~isnan(perfTh) && perfTh > 0
    for ii = 1:length(vectorNames)
        if ~strcmpi(vectorNames{ii},'meanPerfBlockCtrl')
            eval(sprintf('log.%s(log.meanPerfBlockCtrl<perfTh) = [];',vectorNames{ii}));  
        end
    end
end

% apply filters
if isfield(filters,'mazes')
    mazeidx = eval(sprintf('find(log.currMaze%s%d)',filters.mazeLogic,filters.mazes));
else
    mazeidx = 1:length(log.currMaze);
end
if isfield(filters,'power')
    poweridx = eval(sprintf('find(log.power==%1.2f | log.power==0)',filters.power));
else
    poweridx = 1:length(log.currMaze);
end
if iscell(lidx)
    posidx = [];
    for ii = 1:length(log.galvoPosIdxMaster)
        for jj = 1:numel(lidx)
            if isempty(log.galvoPosIdxMaster{ii}) || sum(log.galvoPosIdxMaster{ii}==lidx{jj})==length(lidx{jj})
                posidx = [posidx ii];
            end
        end
    end
    trialIdx = intersect(poweridx,intersect(mazeidx,posidx));
elseif isempty(lidx)
    trialIdx = intersect(poweridx,mazeidx);
else
    posidx = find(log.galvoPosIdxMaster == 0);
    for ii = 1:length(lidx)
        posidx = [posidx find(log.galvoPosIdxMaster == lidx(ii))];
    end
    trialIdx = intersect(poweridx,intersect(mazeidx,posidx));
end

% select relevant trials
for ii = 1:length(vectorNames)
    eval(sprintf('log.%s = log.%s(trialIdx);',vectorNames{ii},vectorNames{ii}));    
end

end