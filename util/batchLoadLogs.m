function logs = batchLoadLogs(fl)

% logs = batchLoadLogs(fl,restrictFlag)
% loads logs in cell array list of file paths "fl" and returns a structure
% array of processed virmen logs


% if files were generated inside spock paths are different
% this is in case i' trying to load locally
if ~isempty(strfind(fl{1},'jukebox'))
  if ~isdir('/jukebox/braininit')
    for iF = 1:numel(fl)
      fl{iF} = ['/Volumes/' fl{iF}(10:end)];
    end
  end
end

ct = 1;
for ii = 1:length(fl)
  load(fl{ii},'logSumm')
  logSumm  = orderfields(logSumm);
  logs(ct) = logSumm;
  ct = ct +1;
end