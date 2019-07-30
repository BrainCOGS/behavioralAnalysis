%PURPLE Purple-white colormap

function map = purple(varargin)
map = [0 0 0 1; 0.25 0 0.25 5; 0.9 0.3 1 4; 0.95 0.85 0.95 0];
map = colormap_helper(map, varargin{:});
