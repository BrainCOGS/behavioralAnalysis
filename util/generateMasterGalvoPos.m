function logs = generateMasterGalvoPos(logs)

% logs = generateMasterGalvoPos(logs)
% add to each log structure in an array a vector containing galvoPosIdx
% referenced to a master grid

for ii = 1:numel(logs)
  grid = logs(ii).info.grid;
  midx = getMasterLocationIdx(grid);
  if iscell(midx)
    logs(ii).galvoPosIdxMaster = cell(size(logs(ii).galvoPosIdx));
    for m = 1:numel(midx)
      thisidx = find(logs(ii).galvoPosIdx==m);
      for jj = 1:numel(thisidx)
        logs(ii).galvoPosIdxMaster{thisidx(jj)} = midx{m};
      end
    end
    for jj = 1:numel(logs(ii).galvoPosIdx)
      if isnan(logs(ii).galvoPosIdx(jj))
        logs(ii).galvoPosIdxMaster{jj} = nan;
      end
    end
  else
    logs(ii).galvoPosIdxMaster = cell(size(logs(ii).galvoPosIdx));
    for m = 1:numel(midx)
      logs(ii).galvoPosIdxMaster(logs(ii).galvoPosIdx==m) = {midx(m)};
    end
    logs(ii).galvoPosIdxMaster(isnan(logs(ii).galvoPosIdx)) = {nan};
  end
  
end