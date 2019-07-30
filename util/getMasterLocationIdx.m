function [idx,mastergrid] = getMasterLocationIdx(coord)

% [idx,mastergrid] = getMasterLocationIdx(coord)
% index of location in a master grid that should span all locations in all
% other grids
% coord is n x 2 matrix with list of [ML AP] coordinates

APs = 4:-.25:-5;
MLs = -4:.25:4;

mastergrid = [];
for ii = 1:numel(APs)
  for jj = 1:numel(MLs)
    mastergrid = [mastergrid; MLs(jj) APs(ii)];
  end
end

if iscell(coord)
  for jj = 1:numel(coord)
    idx{jj} = zeros(size(coord{jj},1),1);
    for ii = 1:numel(idx{jj})
      if coord{jj}(ii,1) == 1.7; coord{jj}(ii,1) = 1.75;
      elseif coord{jj}(ii,1) == -1.7; coord{jj}(ii,1) = -1.75;
      end
      idx{jj}(ii) = find(mastergrid(:,1)==coord{jj}(ii,1) & mastergrid(:,2)==coord{jj}(ii,2));
    end
  end
else
  idx = zeros(size(coord,1),1);
  for ii = 1:numel(idx)
    idx(ii) = find(mastergrid(:,1)==coord(ii,1) & mastergrid(:,2)==coord(ii,2));
  end
end