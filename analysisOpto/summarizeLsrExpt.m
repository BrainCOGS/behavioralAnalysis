function logSumm = summarizeLsrExpt(fp)

% logSumm = summarizeLsrExpt(fp)
% analyzes behavioral data saved in virmen log structure
%
% INPUT
%   fp: cell array with full file paths (typically generated by fileCellArrayLsr.m)
%
% OUTPUT
%   logSumm: data structure with overall performance indicators and trial
%            data
%
% LP feb 16


%% handle files and folders
if iscell(fp.lsr); fp.lsr = fp.lsr{1};       end
if iscell(fp.behav); fp.behav = fp.behav{1}; end

fprintf('   analyzing %s\n',fp.lsr)

% first check if log already exists
bksl = strfind(fp.lsr,'/');
thisname = fp.lsr(bksl(end)+1:end); % file name
% mouse name (try 4 characters first)
temp = regexp(thisname,'[a-z][a-z][0-9][0-9]','match');
if isempty(temp)
    temp = regexp(thisname,'[a-z][a-z][0-9][a-z]','match');
end
if isempty(temp)
    temp = regexp(thisname,'[a-z][a-z][0-9]','match');
end
thismouse = temp{1};
bksl = strfind(thisname,'_');
fn = sprintf('%slogSumm.mat',thisname(1:bksl(end)));

cd(analysisParams.savepath)
if isempty(dir(thismouse))
    mkdir(thismouse)
end
cd(thismouse)
if ~isempty(dir(fn)) % try loading it first
    load(fn,'logSumm')
    if isfield(logSumm,'meanPerfBlockCtrl') && isfield(logSumm,'stimID')
        runFlag = 0;
    else
        runFlag = 1;
        clear logSumm
    end
else
    runFlag = 1;
end

if ~runFlag
    return
else
    % load files from server
    load(fp.behav,'log')
    load(fp.lsr,'lsrlog')
    
    % to more easily access it, concatenate blocks
    lsrtemp = [];
    for ii = 1:length(lsrlog.block)
        lsrtemp = [lsrtemp lsrlog.block(ii).trial(1:end-1)];
    end
    
    %% get general info first
    exptParams = {'mouseID','refIm','currIm','pxlPerMM','refPxl','galvofreq','grid'};
    for ii = 1:length(exptParams)
        eval(sprintf('logSumm.info.%s = lsrlog.block(1).info.exptParams.%s;', exptParams{ii}, exptParams{ii}))
    end
    % deduce genotype from mouse name
    if strcmpi(thismouse(1:2),'vg')
        logSumm.info.genotype = 'vgat';
    else
        logSumm.info.genotype = 'ctrl';
    end
    
    logSumm.info.lsrStartTime = lsrlog.block(1).info.exptParams.startTime;
    logSumm.info.lsrEndTime   = lsrlog.block(1).info.exptParams.endTime;
    logSumm.info.lsrFn        = lsrlog.block(1).info.exptParams.fn;
    logSumm.info.protocol     = log.animal.protocol;
    logSumm.info.laserParams  = log.block(1).laserParams;
    logSumm.info.fn           = fn;
    bksl                      = strfind(fp.behav,'/');
    logSumm.info.behavFn      = fp.behav(bksl(end)+1:end); % file name
    
    temp                      = regexp(fp.lsr,'[0-9]{8,}','match');
    logSumm.info.date         = temp{1}; % file name
    
    mazeInfo = feval(logSumm.info.protocol);
    for ii = 1:length(log.block)
        if ~isempty(log.block(ii).trial)
            logSumm.info.mazes(ii) = log.block(ii).mazeID;
            try
                [~,logSumm.info.mazeLbl{ii}] = formatMazeName(mazeInfo(log.block(ii).mazeID));
            catch
                logSumm.info.mazeLbl{ii} = '';
            end
        end
    end
    logSumm.info.mainMazeID   = log.block(end).mainMazeID;
    
    % first couple of experiments don't have grid label but are all vctx
    % unilateral
    if ~isfield(logSumm.info.laserParams,'gridLabel')
        logSumm.info.laserParams.gridLabel = 'VctxUnilateralGrid.mat';
    end
    
    %% start compiling trials (drop block structure)
    logSumm.keyFrameLabels      = {'cue';'mem';'arm';'rew';'lsrOn';'lsrOff'};
    logSumm.meanPerfBlock       = [];
    logSumm.meanPerfBlockCtrl   = [];
    
    % trial info
    tc = 0; % overall counter to concatenate blocks
    for b = 1:length(log.block) % main maze blocks
        
        % there is a bug that generates empty blocks
        if isempty(log.block(b).trial)
            continue
        end
        
        for t = 1:length(log.block(b).trial) % trials within block
            % which maze
            logSumm.currMaze(tc+t) = log.block(b).mazeID;
            
            % left or right trial
            if strcmpi(char(log.block(b).trial(t).trialType),'L')
                logSumm.trialType(tc+t) = analysisParams.leftCode;
            else
                logSumm.trialType(tc+t) = analysisParams.rightCode;
            end
            
            logSumm.stimID(tc+t) = log.block(b).trial(t).trialID;
            
            % for user-enforced trial ending
            if strcmpi(char(log.block(b).trial(t).choice),'nil') ...
                    && log.block(b).trial(t).time(log.block(b).trial(t).iterations) < 60    
                logSumm = fillwithnan(logSumm,tc,t);
                continue
            end
            
            % laser info
            logSumm.laserON(tc+t)       = log.block(b).trial(t).laserON; % boolean
            if logSumm.laserON(tc+t)
                if isfield(lsrtemp,'power')
                    logSumm.power(tc+t) = lsrtemp(tc+t).power;
                else
                    logSumm.power(tc+t) = logSumm.info.laserParams.power;
                end
            else
                logSumm.power(tc+t) = 0;
            end
            
            thislocation                = unique(log.block(b).trial(t).galvoPosIdx(log.block(b).trial(t).galvoPosIdx>0),'stable'); 
            if isempty(thislocation)
                logSumm.galvoPosIdx(tc+t)   = 0;
            else
                logSumm.galvoPosIdx(tc+t) = thislocation;
            end
            logSumm.laserEpochIdx(tc+t) = log.block(b).trial(t).laserEpoch;
            
            % animal went left or right
            if strcmpi(char(log.block(b).trial(t).choice),'L')
                logSumm.choice(tc+t) = analysisParams.leftCode;
            elseif strcmpi(char(log.block(b).trial(t).choice),'R')
                logSumm.choice(tc+t) = analysisParams.rightCode;
            else
                logSumm.choice(tc+t) = analysisParams.nilCode;
            end
            
            % position, displacement, speed, view angle
            logSumm.excessTravel(tc+t) = log.block(b).trial(t).excessTravel;
            logSumm.pos{tc+t}          = log.block(b).trial(t).position; % position in Maze, cell array; [x y angle]
            logSumm.pos{tc+t}(:,3)     = rad2deg(logSumm.pos{tc+t}(:,3)); % convert view angle to deg
            logSumm.displ{tc+t}        = diff(logSumm.pos{tc+t}); % displacement in Maze, cell array; [x y angle]
            if log.block(b).trial(t).iCueEntry > 0
                stemXYdispl            = logSumm.displ{tc+t}(log.block(b).trial(t).iCueEntry:log.block(b).trial(t).iArmEntry-1,1:2); % total XY displacement in Maze stem
            else
                stemXYdispl            = 0;
            end
            logSumm.stemDispl(tc+t)    = sum(sqrt(sum(stemXYdispl.^2,2))); % total displacement in stem
            if ~strcmpi(char(logSumm.info.protocol),'noworld')
                logSumm.stemDisplNorm(tc+t)= sum(sqrt(sum(stemXYdispl.^2,2)))./...
                    (mazeInfo(logSumm.currMaze(tc+t)).lCue+mazeInfo(logSumm.currMaze(tc+t)).lMemory); % total displacement in stem normalized by stem length
            else
                logSumm.stemDisplNorm(tc+t)= sum(sqrt(sum(stemXYdispl.^2,2)))./...
                    (mazeInfo(1).lCue+mazeInfo(1).lMemory);
            end
            if log.block(b).trial(t).iCueEntry > 0
                if log.block(b).trial(t).iArmEntry > 0
                    logSumm.speedStem(tc+t)    = logSumm.stemDispl(tc+t)/...
                        (log.block(b).trial(t).time(log.block(b).trial(t).iArmEntry)-...
                        log.block(b).trial(t).time(log.block(b).trial(t).iCueEntry)); % speed in maze stem
                else
                    logSumm.speedStem(tc+t)    = logSumm.stemDispl(tc+t)/...
                        (log.block(b).trial(t).time(length(log.block(b).trial(t).position))-...
                        log.block(b).trial(t).time(log.block(b).trial(t).iCueEntry)); % speed in maze stem
                end
            else
                logSumm.speedStem(tc+t)    = nan;
            end

            % cue info
            logSumm.cueOrder{tc+t}  = log.block(b).trial(t).cueCombo; % boolean, [left cues; right cues]
            logSumm.cuePos_R{tc+t}  = log.block(b).trial(t).cuePos{2};
            logSumm.cuePos_L{tc+t}  = log.block(b).trial(t).cuePos{1};
            logSumm.nCues_R(tc+t)   = length(log.block(b).trial(t).cuePos{2});
            logSumm.nCues_L(tc+t)   = length(log.block(b).trial(t).cuePos{1});
            try
                logSumm.cueDur_R{tc+t}  = log.block(b).trial(t).time(log.block(b).trial(t).cueOffset{2}) - ...
                    log.block(b).trial(t).time(log.block(b).trial(t).cueOnset{2});
                logSumm.cueDur_L{tc+t}  = log.block(b).trial(t).time(log.block(b).trial(t).cueOffset{1}) - ...
                    log.block(b).trial(t).time(log.block(b).trial(t).cueOnset{1});
            catch
                logSumm.cueDur_R{tc+t}  = inf; logSumm.cueDur_L{tc+t}  = inf;
            end
            
            % timing info
            logSumm.time{tc+t}          = log.block(b).trial(t).time;
            logSumm.keyFrames{tc+t}     = [log.block(b).trial(t).iCueEntry log.block(b).trial(t).iMemEntry ...
                log.block(b).trial(t).iArmEntry length(log.block(b).trial(t).position) ...
                log.block(b).trial(t).iLaserOn log.block(b).trial(t).iLaserOff]; % corresponds to logSumm.keyFrameLabels; [cue mem arm rew]
            logSumm.trialDur(tc+t)      = logSumm.time{tc+t}(length(log.block(b).trial(t).position));
            if log.block(b).trial(t).iCueEntry > 0
                if log.block(b).trial(t).iMemEntry > 0
                    logSumm.trialDurCue(tc+t)   = logSumm.time{tc+t}(log.block(b).trial(t).iMemEntry) - ...
                        logSumm.time{tc+t}(log.block(b).trial(t).iCueEntry);
                else
                    logSumm.trialDurCue(tc+t)   = logSumm.time{tc+t}(length(log.block(b).trial(t).position)) - ...
                        logSumm.time{tc+t}(log.block(b).trial(t).iCueEntry);
                end
                if log.block(b).trial(t).iMemEntry > 0 && log.block(b).trial(t).iArmEntry > 0
                    logSumm.trialDurMem(tc+t)   = logSumm.time{tc+t}(log.block(b).trial(t).iArmEntry) - ...
                        logSumm.time{tc+t}(log.block(b).trial(t).iMemEntry);
                elseif log.block(b).trial(t).iMemEntry > 0 && log.block(b).trial(t).iArmEntry <= 0
                    logSumm.trialDurMem(tc+t)   = logSumm.time{tc+t}(length(log.block(b).trial(t).position)) - ...
                        logSumm.time{tc+t}(log.block(b).trial(t).iMemEntry);
                else
                    logSumm.trialDurMem(tc+t)   = nan;
                end
            else
                logSumm.trialDurCue(tc+t)   = nan;
                logSumm.trialDurMem(tc+t)   = nan;
            end
            logSumm.trialDurFull(tc+t)  = log.block(b).trial(t).duration;
            try
                fr(tc+t)                = mean(diff(logSumm.time{tc+t}(logSumm.keyFrames{tc+t}(2):logSumm.keyFrames{tc+t}(4))));
            catch
                fr(tc+t)                = nan;
            end
        end
        
        if ~isempty(t);
            % mean performance by block, all trials
            ttype = logSumm.trialType(end-t+1:end);
            ch    = logSumm.choice(end-t+1:end);
            lsr   = logSumm.laserON(end-t+1:end);
            perftemp                = sum(ttype == ch)/t;
            logSumm.meanPerfBlock   = [logSumm.meanPerfBlock ones(1,t)*perftemp];
            
            % mean performance by block, ctrl trials
            perftemp                    = sum(ttype(lsr==0) == ch(lsr==0))/sum(lsr==0);
            logSumm.meanPerfBlockCtrl   = [logSumm.meanPerfBlockCtrl ones(1,t)*perftemp]; 
            
            tc = tc + t; % udpate total number of trials
        end
    end
    
    logSumm.nCues_RminusL = logSumm.nCues_R-logSumm.nCues_L;
    logSumm.ntrials       = length(logSumm.trialType);
    
    % average frame rate during maze navigation (excludes ITI)
    logSumm.frameRate     = nanmean(fr(fr~=0));
%     logSumm = orderfields(logSumm);
    save(sprintf('%s%s/%s',analysisParams.savepath,thismouse,fn),'logSumm')
end

end

%% fill with nan's in case of aborted trial
function logSumm = fillwithnan(logSumm,tc,t)

logSumm.laserON(tc+t)       = nan;
logSumm.power(tc+t)         = nan;
logSumm.galvoPosIdx(tc+t)   = nan;
logSumm.laserEpochIdx(tc+t) = nan;
logSumm.choice(tc+t)        = nan;
logSumm.excessTravel(tc+t)  = nan;
logSumm.pos{tc+t}           = [nan nan nan];
logSumm.displ{tc+t}         = [nan nan nan];
logSumm.stemDispl(tc+t)     = nan;
logSumm.stemDisplNorm(tc+t) = nan;
logSumm.speedStem(tc+t)     = nan;
logSumm.cueOrder{tc+t}      = nan;
logSumm.cuePos_R{tc+t}      = nan;
logSumm.cuePos_L{tc+t}      = nan;
logSumm.nCues_R(tc+t)       = nan;
logSumm.nCues_L(tc+t)       = nan;
logSumm.cueDur_R{tc+t}      = inf;
logSumm.cueDur_L{tc+t}      = inf;
logSumm.time{tc+t}          = nan;
logSumm.keyFrames{tc+t}     = nan;
logSumm.trialDur(tc+t)      = nan;
logSumm.trialDurCue(tc+t)   = nan;
logSumm.trialDurMem(tc+t)   = nan;
logSumm.trialDurFull(tc+t)  = nan;

end

