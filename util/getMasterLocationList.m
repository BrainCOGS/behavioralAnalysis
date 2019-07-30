function masterList = getMasterLocationList

% get list of all unique locations in all grids

temp     = dir(analysisParams.gridpath);
gridList = {temp(3:end).name};

fulllist = {};
for ii = 1:length(gridList)
    temp = locationList(gridList{ii});
    fulllist(end+1:end+length(temp)) = temp;
end

masterList = unique(fulllist);
    