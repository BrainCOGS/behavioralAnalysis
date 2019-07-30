function logs = summarizeLsrExpt_batch(mice,filters,protocol)

% logs = summarizeLsrExpt_batch(mice,filters,protocol)
% applies filters in "filters" data structure to load logs meeting data
% selection criteria. Logs are loaded as a structure array

%% defaults
if nargin < 1; mice      = analysisParams.mice; end
if nargin < 2; filters   = [];                  end
if nargin < 3; protocol  = [];                  end

if ~iscell(mice);    mice     = {mice};         end
if iscell(protocol); protocol = protocol{1};    end

try

if isThisSpock 
  rootdir = analysisParams.pathForSpockBehav; 
  warning('off','all'); 
  addpath(genpath('/usr/people/lpinto/code/tankmousevr/')); 
  warning('on','all'); 
else
  rootdir = analysisParams.serverPathBehav;
end

%% compile list, load logs, update metaLog
if isempty(dir(rootdir))
  fprintf('Server not available. Proceeding with data from local disk\n')
  % go to metalog to filter files
  if nargout > 0
    fl   = filterLogFiles(filters);
    logs = batchLoadLogs(fl);
  end
else
  flp    = fileCellArrayLsr(mice,[],protocol,isThisSpock); % retrieve file list
  
  % create summarized log w/ behav & laser/galvo info
  if numel(flp.lsr) > 0
    if nargout > 0 % concat logs if required, otherwise do one at a time
      logs = {};
      for iFile = 1:numel(flp.lsr)
        fl.lsr        = flp.lsr{iFile};
        fl.behav      = flp.behav{iFile};
        
        [~,thisdate]  = mouseAndDateFromFileName(fl.behav);
        if str2double(thisdate) < 20160328; continue; end % below this date is the 4-m version
        
        try
          logs{end+1} = summarizeVirmenLogLsr(fl);
        catch ME
          displayException(ME)
        end
      end
      
      temp = logs;
      logs = logs{1};
      for iLog = 2:numel(temp); logs(iLog) = temp{iLog}; end
      
      mkLsrMetaLog(logs); % update metalog file;
      
    else
      for iFile = 1:numel(flp.lsr)
        fl.lsr        = flp.lsr{iFile};
        fl.behav      = flp.behav{iFile};
        
        [~,thisdate]  = mouseAndDateFromFileName(fl.behav);
        if str2double(thisdate) < 20160328; continue; end % below this date is the 4-m version
        
        try
          mkLsrMetaLog(summarizeVirmenLogLsr(fl)); % update metalog file;
        catch ME
          displayException(ME)
        end
      end
      
    end
  else
    logs = [];
  end
  
end
catch ME
  displayException(ME)
end
