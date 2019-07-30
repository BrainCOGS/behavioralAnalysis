%RED Black-red-orange colormap

function map = red(varargin)
map = [0 0 0 1; 0.25 0 0 4; 0.8 0 0 2; 1 0 0 4; 0.9 0.5 0 2; 0.9 0.7 0 5];
map = colormap_helper(map, varargin{:});
